#include "gaussian_smc_adaptor.h"
#include "proposal_kernel.h"
#include "smc_output.h"
#include "particles.h"
#include "move_output.h"
#include "ensemble_kalman_output.h"
#include "ensemble.h"

namespace ilike
{

GaussianSMCAdaptor::GaussianSMCAdaptor()
:SMCAdaptor()
{
  this->constant = -1.0;
}

GaussianSMCAdaptor::GaussianSMCAdaptor(const std::string &variable_in,
                                       double constant_in)
:SMCAdaptor()
{
  this->variable = variable_in;
  this->constant = constant_in;
}

GaussianSMCAdaptor::~GaussianSMCAdaptor()
{
}

GaussianSMCAdaptor::GaussianSMCAdaptor(const GaussianSMCAdaptor &another)
:SMCAdaptor(another)
{
  this->make_copy(another);
}

void GaussianSMCAdaptor::operator=(const GaussianSMCAdaptor &another)
{
  if(this == &another)
    return;
  SMCAdaptor::operator=(another);
  this->make_copy(another);
}

void GaussianSMCAdaptor::make_copy(const GaussianSMCAdaptor &another)
{
  this->variable = another.variable;
  this->constant = another.constant;
}

SMCAdaptor* GaussianSMCAdaptor::duplicate() const
{
  return( new GaussianSMCAdaptor(*this));
}

void GaussianSMCAdaptor::adapt_from(const arma::mat &points,
                                    const arma::colvec &normalised_weights)
{
  if (this->proposal==NULL)
    Rcpp::stop("SMC adaptor: not attached to a proposal.");
  if (this->proposal->get_transform()!=NULL)
    Rcpp::stop("SMC adaptor: not available for a proposal on a transformed space.");
  GaussianProposalInfo* info = this->proposal->get_gaussian_proposal_info(this->variable);
  if (info==NULL)
    Rcpp::stop("SMC adaptor: the proposal has no Gaussian covariance for variable " + this->variable + " to adapt.");
  PrincipalComponentsInfo* components = this->proposal->get_principal_components(this->variable);
  arma::colvec* component_scales = this->proposal->get_component_scales(this->variable);

  size_t d = points.n_cols;
  if (points.n_rows<2)
    return;
  if (info->get_covariance().n_rows!=d)
    Rcpp::stop("SMC adaptor: the proposal's covariance for " + this->variable + " does not match the variable's dimension.");

  // weighted mean and covariance (the unbiased form for normalised weights)
  arma::rowvec mean = normalised_weights.t()*points;
  arma::mat centred = points.each_row() - mean;
  double denominator = 1.0 - arma::accu(arma::square(normalised_weights));
  if (denominator<=0.0)
    return; // a single particle carries all the weight
  arma::mat covariance = (centred.t()*(centred.each_col() % normalised_weights))/denominator;
  covariance = (covariance + covariance.t())/2.0;

  if (this->proposal->proposal_uses_gaussian_mean())
    info->get_mean() = mean.t();

  if (covariance.is_sympd())
  {
    info->set_covariance(covariance);
    if (components!=NULL)
    {
      arma::colvec values;
      arma::mat vectors;
      if (arma::eig_sym(values, vectors, covariance))
      {
        size_t m = components->directions.n_cols;
        components->directions = arma::fliplr(vectors.tail_cols(m));
        components->eigenvalues = arma::flipud(values.tail(m));
      }
    }
  }

  if (this->constant>0.0)
  {
    double c2 = this->constant*this->constant;
    if (components!=NULL)
      components->scales.fill(c2);
    else if (component_scales!=NULL)
      component_scales->fill(c2);
    else
      info->get_double_scale() = c2/double(d);
  }
}

void GaussianSMCAdaptor::smc_adapt(SMCOutput* current_state)
{
  Particles &generation = current_state->back();
  arma::mat points;
  if (generation.is_flat())
  {
    std::vector<std::string> variables(1, this->variable);
    points = generation.flat_points.matrix_of(variables);
  }
  else
  {
    size_t n = generation.particles.size();
    for (size_t i=0; i<n; ++i)
    {
      arma::colvec point = generation.particles[i]->back().parameters.get_colvec(this->variable);
      if (i==0)
        points.set_size(n, point.n_elem);
      points.row(i) = point.t();
    }
  }
  if (points.n_rows==0)
    return;
  // the weights for the current target, before resampling
  arma::colvec log_weights = generation.unnormalised_log_weights;
  arma::colvec weights;
  if (log_weights.n_elem==points.n_rows)
  {
    weights = arma::exp(log_weights - log_weights.max());
    if (!weights.is_finite() || (arma::accu(weights)<=0.0))
      return;
    weights = weights/arma::accu(weights);
  }
  else
  {
    weights = arma::colvec(points.n_rows, arma::fill::value(1.0/double(points.n_rows)));
  }
  this->adapt_from(points, weights);
}

void GaussianSMCAdaptor::ensemble_adapt(EnsembleKalmanOutput* current_state)
{
  // ensemble members are equally weighted
  Ensemble &ensemble = current_state->back();
  size_t n = ensemble.members.size();
  if (n==0)
    return;
  arma::mat points;
  for (size_t i=0; i<n; ++i)
  {
    arma::colvec point = ensemble.members[i]->back().parameters.get_colvec(this->variable);
    if (i==0)
      points.set_size(n, point.n_elem);
    points.row(i) = point.t();
  }
  this->adapt_from(points, arma::colvec(n, arma::fill::value(1.0/double(n))));
}

}
