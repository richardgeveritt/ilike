#include "gaussian_mcmc_adaptor.h"
#include "proposal_kernel.h"
#include "particle.h"

namespace ilike
{

double mirror_default_target_acceptance(size_t d)
{
  // optimal acceptance rates found by experiment (dev/plan.md 3.2, C171)
  const double dimensions[] = {1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0, 100.0};
  const double rates[] = {0.50, 0.42, 0.36, 0.33, 0.26, 0.25, 0.23, 0.23};
  const size_t n = 8;
  double x = double(d);
  if (x<=dimensions[0])
    return rates[0];
  if (x>=dimensions[n-1])
    return rates[n-1];
  size_t k = 0;
  while (dimensions[k+1]<x)
    k = k + 1;
  double w = (log(x) - log(dimensions[k]))/(log(dimensions[k+1]) - log(dimensions[k]));
  return (1.0-w)*rates[k] + w*rates[k+1];
}

GaussianMCMCAdaptor::GaussianMCMCAdaptor()
:MCMCAdaptor()
{
  this->method = GaussianAdaptationMethod::SCALED_AM;
  this->target_acceptance = 0.234;
  this->gain_exponent = 0.6;
  this->gain_constant = 1.0;
  this->refresh_every = 0;
  this->number_of_updates = 0;
  this->writes_mean = false;
  this->averaging = false;
  this->bound_to = NULL;
}

GaussianMCMCAdaptor::GaussianMCMCAdaptor(const std::vector<std::string> &variable_names_in,
                                         GaussianAdaptationMethod method_in,
                                         double target_acceptance_in,
                                         double gain_exponent_in,
                                         double gain_constant_in,
                                         size_t refresh_every_in,
                                         bool averaging_in)
:MCMCAdaptor()
{
  if ( (target_acceptance_in<0.0) || (target_acceptance_in>=1.0) )
    Rcpp::stop("MCMC adaptor: the target acceptance rate must be between 0 and 1.");
  if ( (gain_exponent_in<=0.5) || (gain_exponent_in>1.0) )
    Rcpp::stop("MCMC adaptor: the gain exponent must be in (0.5, 1].");
  if (gain_constant_in<=0.0)
    Rcpp::stop("MCMC adaptor: the gain constant must be positive.");
  for (auto i=variable_names_in.begin(); i!=variable_names_in.end(); ++i)
  {
    VariableState state;
    state.variable = *i;
    state.info = NULL;
    state.component_scales = NULL;
    state.components = NULL;
    state.log_scale = 0.0;
    state.target = target_acceptance_in;
    state.started = false;
    this->states.push_back(state);
  }
  this->method = method_in;
  this->target_acceptance = target_acceptance_in;
  this->gain_exponent = gain_exponent_in;
  this->gain_constant = gain_constant_in;
  this->refresh_every = refresh_every_in;
  this->number_of_updates = 0;
  this->writes_mean = false;
  this->averaging = averaging_in;
  this->bound_to = NULL;
}

GaussianMCMCAdaptor::~GaussianMCMCAdaptor()
{
}

GaussianMCMCAdaptor::GaussianMCMCAdaptor(const GaussianMCMCAdaptor &another)
:MCMCAdaptor(another)
{
  this->make_copy(another);
}

void GaussianMCMCAdaptor::operator=(const GaussianMCMCAdaptor &another)
{
  if(this == &another)
    return;
  MCMCAdaptor::operator=(another);
  this->make_copy(another);
}

void GaussianMCMCAdaptor::make_copy(const GaussianMCMCAdaptor &another)
{
  this->states = another.states;
  this->method = another.method;
  this->target_acceptance = another.target_acceptance;
  this->gain_exponent = another.gain_exponent;
  this->gain_constant = another.gain_constant;
  this->refresh_every = another.refresh_every;
  this->number_of_updates = another.number_of_updates;
  this->writes_mean = another.writes_mean;
  this->averaging = another.averaging;
  // the pointers are into the original kernel: looked up again in the copy's
  // kernel before the next update
  this->bound_to = NULL;
  for (auto i=this->states.begin(); i!=this->states.end(); ++i)
  {
    i->info = NULL;
    i->component_scales = NULL;
    i->components = NULL;
  }
}

MCMCAdaptor* GaussianMCMCAdaptor::duplicate() const
{
  return( new GaussianMCMCAdaptor(*this));
}

bool GaussianMCMCAdaptor::uses_global_scale() const
{
  return (this->method==GaussianAdaptationMethod::AM) ||
    (this->method==GaussianAdaptationMethod::RB_AM) ||
    (this->method==GaussianAdaptationMethod::SCALED_AM) ||
    (this->method==GaussianAdaptationMethod::SCALE);
}

bool GaussianMCMCAdaptor::estimates_moments() const
{
  return this->method!=GaussianAdaptationMethod::SCALE;
}

void GaussianMCMCAdaptor::bind()
{
  if (this->proposal==NULL)
    Rcpp::stop("MCMC adaptor: not attached to a proposal.");
  this->writes_mean = this->proposal->proposal_uses_gaussian_mean();
  for (auto i=this->states.begin(); i!=this->states.end(); ++i)
  {
    i->info = this->proposal->get_gaussian_proposal_info(i->variable);
    if (i->info==NULL)
      Rcpp::stop("MCMC adaptor: the proposal has no Gaussian covariance for variable " + i->variable + " to adapt.");
    i->component_scales = this->proposal->get_component_scales(i->variable);
    if ( (this->method==GaussianAdaptationMethod::COMPONENTWISE_AM) && (i->component_scales==NULL) )
      Rcpp::stop("MCMC adaptor: ilike::componentwise_am needs the proposal ilike::componentwise_rw.");
    i->components = this->proposal->get_principal_components(i->variable);
    if ( ((this->method==GaussianAdaptationMethod::PC_EIGEN) || (this->method==GaussianAdaptationMethod::PC_ONLINE)) && (i->components==NULL) )
      Rcpp::stop("MCMC adaptor: ilike::pc_am and ilike::online_pc_am need the proposal ilike::pc_rw.");
  }
  this->bound_to = this->proposal;
}

void GaussianMCMCAdaptor::refresh_principal_components(VariableState &s,
                                                       const arma::mat &covariance)
{
  arma::mat symmetric = (covariance + covariance.t())/2.0;
  size_t m = s.components->directions.n_cols;
  if (this->method==GaussianAdaptationMethod::PC_EIGEN)
  {
    arma::colvec values;
    arma::mat vectors;
    if (!arma::eig_sym(values, vectors, symmetric))
      return;
    // ascending: the first m principal directions are the last m
    s.components->directions = arma::fliplr(vectors.tail_cols(m));
    s.components->eigenvalues = arma::flipud(values.tail(m));
  }
  else
  {
    // one step of the online PCA recursion (Section 5.3.2; the constrained
    // EM of Ahn & Oh 2003 with sigma^2 = 0), with L(.) and U(.) the lower and
    // upper triangular parts:
    //   W <- Sigma W L(W'W)^{-T} [ U( L(W'W)^{-1} W' Sigma W L(W'W)^{-T} ) ]^{-1}
    // Its fixed points are W = (eigenvectors) D, D diagonal: the columns
    // converge to the principal directions in order.
    const arma::mat &W = s.pca_loadings;
    arma::mat lower = arma::trimatl(W.t()*W);
    arma::mat lower_inverse;
    if (!arma::inv(lower_inverse, arma::trimatl(lower)))
      return;
    arma::mat projected = symmetric*W*lower_inverse.t();
    arma::mat inner = lower_inverse*W.t()*projected;
    arma::mat upper_inverse;
    if (!arma::inv(upper_inverse, arma::trimatu(inner)))
      return;
    arma::mat updated = projected*upper_inverse;
    if (!updated.is_finite())
      return;
    s.pca_loadings = updated;
    // the kernel's directions: orthonormal (so that the direction moved can
    // be found from Y - X), in the order of the columns of W
    arma::mat Q, R;
    if (!arma::qr_econ(Q, R, updated))
      return;
    for (size_t l=0; l<m; ++l)
    {
      if (R(l,l)<0.0)
        Q.col(l) = -Q.col(l);
    }
    s.components->directions = Q;
    for (size_t l=0; l<m; ++l)
      s.components->eigenvalues[l] = arma::as_scalar(Q.col(l).t()*symmetric*Q.col(l));
  }
  for (size_t l=0; l<m; ++l)
  {
    if (s.components->eigenvalues[l]<=0.0)
      s.components->eigenvalues[l] = 1e-12;
  }
}

void GaussianMCMCAdaptor::adapt(const AdaptationStep &step)
{
  if (this->bound_to!=this->proposal)
    this->bind();

  this->number_of_updates = this->number_of_updates + 1;
  double gain = this->gain_constant*pow(double(this->number_of_updates+1), -this->gain_exponent);
  double alpha = step.acceptance_probability;

  for (auto s=this->states.begin(); s!=this->states.end(); ++s)
  {
    // the kernel proposes on its own (possibly transformed) space
    arma::colvec previous = step.previous->get_transformed_parameters(this->proposal).get_colvec(s->variable);

    if (!s->started)
    {
      // mu_0 = X_0 (for a mirror proposal, its mean), Sigma_0 = the
      // proposal's covariance, lambda_0 = its scale (2.38^2/d for AM and
      // RB_AM on a random walk, whose scale is not adapted)
      s->mean = previous;
      if (this->writes_mean && (s->info->get_mean().n_elem==previous.n_elem))
        s->mean = s->info->get_mean();
      s->covariance = s->info->get_covariance();
      if (s->covariance.n_rows!=previous.n_elem)
        Rcpp::stop("MCMC adaptor: the proposal's covariance for " + s->variable + " does not match the variable's dimension.");
      if ( ((this->method==GaussianAdaptationMethod::AM) || (this->method==GaussianAdaptationMethod::RB_AM)) && !this->writes_mean )
      {
        s->info->get_double_scale() = 2.38*2.38/double(previous.n_elem);
      }
      s->log_scale = log(s->info->get_double_scale());
      if (s->target==0.0)
        s->target = mirror_default_target_acceptance(previous.n_elem);
      if (s->component_scales!=NULL)
        s->log_component_scales = arma::log(*s->component_scales);
      if (s->components!=NULL)
      {
        s->log_component_scales = arma::log(s->components->scales);
        // W starts at the kernel's directions scaled by sqrt(rho_l)
        s->pca_loadings = s->components->directions*arma::diagmat(arma::sqrt(s->components->eigenvalues));
      }
      // theta-bar_0 = theta_0
      s->average_log_scale = s->log_scale;
      s->average_mean = s->mean;
      s->average_covariance = s->covariance;
      s->average_log_component_scales = s->log_component_scales;
      s->started = true;
    }

    if ( (this->method==GaussianAdaptationMethod::SCALED_AM) || (this->method==GaussianAdaptationMethod::SCALE) )
    {
      s->log_scale = s->log_scale + gain*(alpha - s->target);
    }

    if (this->method==GaussianAdaptationMethod::COMPONENTWISE_AM)
    {
      // (32): the component moved is where the proposal differs from X_i (a
      // zero increment has probability zero; then no scale is updated)
      arma::colvec increment = step.proposed->get_transformed_parameters(this->proposal).get_colvec(s->variable) - previous;
      arma::uvec moved = arma::find(increment!=0.0, 1);
      if (moved.n_elem==1)
      {
        size_t k = moved[0];
        s->log_component_scales[k] = s->log_component_scales[k] + gain*(alpha - s->target);
      }
    }

    if ( (this->method==GaussianAdaptationMethod::PC_EIGEN) || (this->method==GaussianAdaptationMethod::PC_ONLINE) )
    {
      // Algorithm 8: the direction moved is the one the increment lies along
      // (the kernel's directions are orthonormal)
      arma::colvec increment = step.proposed->get_transformed_parameters(this->proposal).get_colvec(s->variable) - previous;
      if (arma::norm(increment)>0.0)
      {
        arma::colvec projections = arma::abs(s->components->directions.t()*increment);
        size_t l = projections.index_max();
        s->log_component_scales[l] = s->log_component_scales[l] + gain*(alpha - s->target);
      }
    }

    if ( (this->method==GaussianAdaptationMethod::AM) || (this->method==GaussianAdaptationMethod::SCALED_AM) ||
         (this->method==GaussianAdaptationMethod::COMPONENTWISE_AM) ||
         (this->method==GaussianAdaptationMethod::PC_EIGEN) || (this->method==GaussianAdaptationMethod::PC_ONLINE) )
    {
      // (27)/(31): with the mean before this update
      arma::colvec difference = step.current->get_transformed_parameters(this->proposal).get_colvec(s->variable) - s->mean;
      s->mean = s->mean + gain*difference;
      s->covariance = s->covariance + gain*(difference*difference.t() - s->covariance);
    }
    else if (this->method==GaussianAdaptationMethod::RB_AM)
    {
      // (30): the expectation of the update over the accept/reject step
      arma::colvec to_proposed = step.proposed->get_transformed_parameters(this->proposal).get_colvec(s->variable) - s->mean;
      arma::colvec to_previous = previous - s->mean;
      s->mean = s->mean + gain*(alpha*to_proposed + (1.0-alpha)*to_previous);
      s->covariance = s->covariance + gain*(alpha*to_proposed*to_proposed.t() + (1.0-alpha)*to_previous*to_previous.t() - s->covariance);
    }

    // with averaging: theta-bar_n = theta-bar_{n-1} + (theta_n - theta-bar_{n-1})/(n+1)
    if (this->averaging)
    {
      double weight = 1.0/double(this->number_of_updates + 1);
      s->average_log_scale = s->average_log_scale + weight*(s->log_scale - s->average_log_scale);
      if (s->log_component_scales.n_elem>0)
        s->average_log_component_scales = s->average_log_component_scales + weight*(s->log_component_scales - s->average_log_component_scales);
      if (this->estimates_moments())
      {
        s->average_mean = s->average_mean + weight*(s->mean - s->average_mean);
        s->average_covariance = s->average_covariance + weight*(s->covariance - s->average_covariance);
      }
    }

    // what the kernel uses: the iterates, or their averages
    double used_log_scale = this->averaging ? s->average_log_scale : s->log_scale;
    const arma::colvec &used_log_component_scales = this->averaging ? s->average_log_component_scales : s->log_component_scales;
    const arma::colvec &used_mean = this->averaging ? s->average_mean : s->mean;
    const arma::mat &used_covariance = this->averaging ? s->average_covariance : s->covariance;

    if ( (this->method==GaussianAdaptationMethod::SCALED_AM) || (this->method==GaussianAdaptationMethod::SCALE) )
      s->info->get_double_scale() = exp(used_log_scale);
    if (s->component_scales!=NULL)
      *s->component_scales = arma::exp(used_log_component_scales);
    if (s->components!=NULL)
      s->components->scales = arma::exp(used_log_component_scales);

    if (this->writes_mean && this->estimates_moments())
      s->info->get_mean() = used_mean;

    if (this->estimates_moments())
    {
      size_t every = this->refresh_every;
      if (every==0)
        every = (this->method==GaussianAdaptationMethod::PC_ONLINE) ? 1 : s->mean.n_elem;
      if ((this->number_of_updates % every)==0)
      {
        if (s->components!=NULL)
        {
          this->refresh_principal_components(*s, used_covariance);
        }
        else
        {
          // each update is a convex combination of a positive definite
          // matrix and a positive semi-definite one (gain < 1), so Sigma
          // stays positive definite in exact arithmetic; if rounding breaks
          // that, the proposal keeps its previous covariance until the next
          // refresh
          arma::mat symmetric = (used_covariance + used_covariance.t())/2.0;
          if (symmetric.is_sympd())
            s->info->set_covariance(symmetric);
        }
      }
    }
  }
}

void GaussianMCMCAdaptor::append_state(std::vector<std::string> *names,
                                       std::vector<double> *values,
                                       bool final,
                                       const std::string &prefix) const
{
  for (auto s=this->states.begin(); s!=this->states.end(); ++s)
  {
    if (!s->started)
      continue;
    std::string base = prefix + s->variable + ".";
    size_t d = s->mean.n_elem;
    // the values the kernel uses (the averages, with averaging)
    double log_scale = this->averaging ? s->average_log_scale : s->log_scale;
    const arma::colvec &log_component_scales = this->averaging ? s->average_log_component_scales : s->log_component_scales;
    const arma::colvec &mean = this->averaging ? s->average_mean : s->mean;
    const arma::mat &covariance = this->averaging ? s->average_covariance : s->covariance;
    if (this->uses_global_scale())
    {
      if (names!=NULL)
        names->push_back(base + "scale");
      if (values!=NULL)
        values->push_back(exp(log_scale));
    }
    if (s->component_scales!=NULL)
    {
      for (size_t k=0; k<log_component_scales.n_elem; ++k)
      {
        if (names!=NULL)
          names->push_back(base + "component_scale[" + std::to_string(k+1) + "]");
        if (values!=NULL)
          values->push_back(exp(log_component_scales[k]));
      }
    }
    if (s->components!=NULL)
    {
      for (size_t l=0; l<log_component_scales.n_elem; ++l)
      {
        if (names!=NULL)
          names->push_back(base + "direction_scale[" + std::to_string(l+1) + "]");
        if (values!=NULL)
          values->push_back(exp(log_component_scales[l]));
      }
    }
    if (this->estimates_moments())
    {
      for (size_t k=0; k<d; ++k)
      {
        if (names!=NULL)
          names->push_back(base + "mean[" + std::to_string(k+1) + "]");
        if (values!=NULL)
          values->push_back(mean[k]);
      }
    }
    if (!final)
      continue;
    if (this->estimates_moments())
    {
      // the upper triangle, row by row
      for (size_t i=0; i<d; ++i)
      {
        for (size_t j=i; j<d; ++j)
        {
          if (names!=NULL)
            names->push_back(base + "covariance[" + std::to_string(i+1) + "," + std::to_string(j+1) + "]");
          if (values!=NULL)
            values->push_back(covariance(i,j));
        }
      }
    }
    if (s->components!=NULL)
    {
      size_t m = s->components->directions.n_cols;
      for (size_t l=0; l<m; ++l)
      {
        if (names!=NULL)
          names->push_back(base + "eigenvalue[" + std::to_string(l+1) + "]");
        if (values!=NULL)
          values->push_back(s->components->eigenvalues[l]);
      }
      for (size_t l=0; l<m; ++l)
      {
        for (size_t i=0; i<d; ++i)
        {
          if (names!=NULL)
            names->push_back(base + "direction[" + std::to_string(i+1) + "," + std::to_string(l+1) + "]");
          if (values!=NULL)
            values->push_back(s->components->directions(i,l));
        }
      }
    }
  }
}

}
