#include "generic_measurement_covariance_estimator.h"
#include "generic_measurement_covariance_estimator_output.h"

namespace ilike
{
GenericMeasurementCovarianceEstimator::GenericMeasurementCovarianceEstimator()
:MeasurementCovarianceEstimator()
{
}

GenericMeasurementCovarianceEstimator::GenericMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                                             size_t* seed_in,
                                                                             Data* data_in,
                                                                             std::shared_ptr<Transform> transform_in,
                                                                             std::shared_ptr<Transform> summary_statistics_in,
                                                                             SimulateModelPtr simulator_in,
                                                                             const std::vector<std::string> &measurement_variables_in)
: MeasurementCovarianceEstimator(rng_in,
                                 seed_in,
                                 data_in,
                                 transform_in,
                                 summary_statistics_in,
                                 measurement_variables_in)
{
  this->simulator = simulator_in;
}

GenericMeasurementCovarianceEstimator::~GenericMeasurementCovarianceEstimator()
{
}

GenericMeasurementCovarianceEstimator::GenericMeasurementCovarianceEstimator(const GenericMeasurementCovarianceEstimator &another)
:MeasurementCovarianceEstimator(another)
{
  this->make_copy(another);
}

void GenericMeasurementCovarianceEstimator::operator=(const GenericMeasurementCovarianceEstimator &another)
{
  if(this == &another)
    return;
  
  MeasurementCovarianceEstimator::operator=(another);
  this->make_copy(another);
}

void GenericMeasurementCovarianceEstimator::make_copy(const GenericMeasurementCovarianceEstimator &another)
{
  this->dimension = another.dimension;
  this->simulator = another.simulator;
  this->Cygivenx = another.Cygivenx;
  this->gaussian_simulator = another.gaussian_simulator;
}

MeasurementCovarianceEstimator* GenericMeasurementCovarianceEstimator::duplicate() const
{
  return( new GenericMeasurementCovarianceEstimator(*this));
}

MeasurementCovarianceEstimatorOutput* GenericMeasurementCovarianceEstimator::initialise_measurement_covariance_estimator()
{
  MeasurementCovarianceEstimatorOutput* output = new GenericMeasurementCovarianceEstimatorOutput(this);
  return output;
}

MeasurementCovarianceEstimatorOutput* GenericMeasurementCovarianceEstimator::initialise_measurement_covariance_estimator(const Parameters &conditioned_on_parameters)
{
  MeasurementCovarianceEstimatorOutput* output = new GenericMeasurementCovarianceEstimatorOutput(this);
  return output;
}

void GenericMeasurementCovarianceEstimator::setup()
{
  this->setup_measurement_variables();
}

void GenericMeasurementCovarianceEstimator::setup(const Parameters &parameters)
{
  this->setup_measurement_variables(parameters);
}

void GenericMeasurementCovarianceEstimator::setup_measurement_variables()
{
  Data dummy_data = this->simulator(*this->rng,Parameters());
  std::vector<arma::colvec> means;
  means.reserve(this->measurement_variables.size());
  std::vector<arma::mat> covs;
  covs.reserve(this->measurement_variables.size());
  for (size_t i=0; i<this->measurement_variables.size(); ++i)
  {
    arma::colvec current_vector = arma::vectorise(dummy_data[this->measurement_variables[i]]);
    means.push_back(arma::zeros<arma::colvec>(current_vector.n_rows));
    // N(0, I), as in setup_measurement_variables(parameters) (a zero matrix
    // here is not positive definite, B20).
    covs.push_back(arma::eye(current_vector.n_rows,current_vector.n_rows));
  }
  this->gaussian_simulator = GaussianIndependentProposalKernel(this->measurement_variables,
                                                               means,
                                                               covs);
}

void GenericMeasurementCovarianceEstimator::setup_measurement_variables(const Parameters &conditioned_on_parameters)
{
  Data dummy_data = this->simulator(*this->rng,conditioned_on_parameters);
  std::vector<arma::colvec> means;
  means.reserve(this->measurement_variables.size());
  std::vector<arma::mat> covs;
  covs.reserve(this->measurement_variables.size());
  for (size_t i=0; i<this->measurement_variables.size(); ++i)
  {
    arma::colvec current_vector = arma::vectorise(dummy_data[this->measurement_variables[i]]);
    means.push_back(arma::zeros<arma::colvec>(current_vector.n_rows));
    covs.push_back(arma::eye(current_vector.n_rows,current_vector.n_rows));
  }
  this->gaussian_simulator = GaussianIndependentProposalKernel(this->measurement_variables,
                                                               means,
                                                               covs);
}

arma::mat GenericMeasurementCovarianceEstimator::get_Cygivenx() const
{
  return this->Cygivenx;
}

arma::mat GenericMeasurementCovarianceEstimator::get_adjustment(const arma::mat &Zf,
                                                                const arma::mat &Dhathalf,
                                                                const arma::mat &P,
                                                                const arma::mat &Vtranspose,
                                                                const arma::mat &Yhat,
                                                                double inverse_incremental_temperature) const
{
  // follows https://arxiv.org/abs/2006.02941
  arma::mat I;
  I.eye(Vtranspose.n_cols,Vtranspose.n_cols);
  
  arma::mat for_eig = Vtranspose*(arma::inv_sympd(I + Yhat*arma::inv_sympd((inverse_incremental_temperature-1.0)*this->get_Cygivenx())*Yhat.t()))*Vtranspose.t();
  
  for_eig = (for_eig+for_eig.t())/2.0;
  
  arma::mat U;
  arma::vec diagD;
  arma::eig_sym(diagD,U,for_eig);
  arma::mat Dsqrt(diagD.n_elem,diagD.n_elem);
  Dsqrt.diag() = arma::sqrt(diagD);
  
  return P*Dhathalf*U*Dsqrt*arma::pinv(Dhathalf)*P.t();
}

arma::mat GenericMeasurementCovarianceEstimator::get_sqrt_adjustment(const arma::mat &Cxy,
                                                                     const arma::mat &Cyy,
                                                                     double inverse_incremental_temperature) const
{
  Rcpp::stop("GenericMeasurementCovarianceEstimator::get_sqrt_adjustment - not yet implemented.");
  
}

bool GenericMeasurementCovarianceEstimator::need_Cxx() const
{
  return true;
}

void GenericMeasurementCovarianceEstimator::find_Cygivenx(const arma::mat &inv_Cxx,
                                                          const arma::mat &Cxy,
                                                          const arma::mat &Cyy)
{
  this->Cygivenx = Cyy - Cxy.t()*inv_Cxx*Cxy;
}

arma::mat GenericMeasurementCovarianceEstimator::get_unconditional_measurement_covariance(const arma::mat &Cyy,
                                                                                          double inverse_incremental_temperature) const
{
  return Cyy + (inverse_incremental_temperature-1.0)*this->Cygivenx;
}

Parameters GenericMeasurementCovarianceEstimator::simulate(const Parameters &current_state)
{
  return this->simulator(*this->rng,
                         current_state);
}

arma::colvec GenericMeasurementCovarianceEstimator::gaussian_simulate()
{
  return this->gaussian_simulator.independent_simulate(*this->rng).get_colvec(this->measurement_variables);
}

void GenericMeasurementCovarianceEstimator::change_data()
{
  this->current_data = this->data;
  if (this->measurement_variables.size()>0)
  {
    this->measurement = (*this->current_data)[this->measurement_variables[0]].as_col();
    
    for (size_t i=1; i<this->measurement_variables.size(); ++i)
    {
      this->measurement = join_cols(this->measurement,(*this->current_data)[this->measurement_variables[i]].as_col());
    }
  }
}

void GenericMeasurementCovarianceEstimator::change_data(std::shared_ptr<Data> new_data)
{
  this->current_data = new_data.get();
  if (this->measurement_variables.size()>0)
  {
    this->measurement = (*this->current_data)[this->measurement_variables[0]].as_col();
    
    if (this->measurement_variables.size()>0)
    {
      for (size_t i=1; i<this->measurement_variables.size(); ++i)
      {
        this->measurement = join_cols(this->measurement,(*this->current_data)[this->measurement_variables[i]].as_col());
      }
    }
  }
}

void GenericMeasurementCovarianceEstimator::precompute_gaussian_covariance(double inverse_incremental_temperature,
                                                                           arma::mat &inv_sigma_precomp,
                                                                           double &log_det_precomp)
{
  inv_sigma_precomp = arma::inv_sympd(inverse_incremental_temperature*this->Cygivenx);
  log_det_precomp = arma::log_det_sympd(inverse_incremental_temperature*this->Cygivenx);
}

}
