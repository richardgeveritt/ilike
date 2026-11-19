#include "gaussian_smc_adaptor.h"
#include "proposal_kernel.h"
#include "smc_output.h"
#include "particles.h"
#include "move_output.h"
#include "ensemble_kalman_output.h"
#include "ensemble.h"
#include "covariance_tools.h"

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
                                    const arma::colvec &normalised_weights,
                                    size_t d)
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

  if (points.n_rows<2)
    return;
  // a 1 x 1 covariance c on a block of several elements is c I (C205)
  if ((info->get_covariance().n_rows==1) && (d>1))
    info->set_covariance(info->get_covariance()(0,0)*arma::eye(d, d));
  check_proposal_size(this->variable, d, info->get_covariance().n_rows, "covariance");

  // the moved variable (the first d columns) given all the other parameters:
  // its weighted mean and conditional covariance (shrunk if singular)
  arma::rowvec mean;
  arma::mat covariance;
  bool estimated = conditional_moments(points, normalised_weights, d, mean, covariance);
  if (!estimated)
  {
    // (the mean alone, when the covariance cannot be estimated)
    arma::mat full_covariance;
    if (!weighted_covariance(points.cols(0, d - 1), normalised_weights, mean, full_covariance))
      return;
  }

  if (this->proposal->proposal_uses_gaussian_mean())
    info->get_mean() = mean.t();

  if (estimated && covariance.is_sympd())
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
  // the moved block's variables first (several for a stacked variable)
  std::vector<std::string> names = split_stacked_name(this->variable);
  std::vector<size_t> offsets;
  arma::mat points = population_points(generation, names, offsets);
  if (points.n_rows==0)
    return;
  arma::colvec weights;
  if (!population_weights(generation, points.n_rows, weights))
    return;
  size_t K = split_stacked_name(this->variable).size();
  size_t d = (offsets.size() > K) ? offsets[K] : points.n_cols;
  this->adapt_from(points, weights, d);
}

void GaussianSMCAdaptor::ensemble_adapt(EnsembleKalmanOutput* current_state)
{
  // ensemble members are equally weighted
  Ensemble &ensemble = current_state->back();
  // the moved block's variables first (several for a stacked variable)
  std::vector<std::string> names = split_stacked_name(this->variable);
  std::vector<size_t> offsets;
  arma::mat points = ensemble_points(ensemble, names, offsets);
  if (points.n_rows==0)
    return;
  size_t K = split_stacked_name(this->variable).size();
  size_t d = (offsets.size() > K) ? offsets[K] : points.n_cols;
  this->adapt_from(points, arma::colvec(points.n_rows, arma::fill::value(1.0/double(points.n_rows))), d);
}

}
