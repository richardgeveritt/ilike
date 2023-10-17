#include "mixed_generic_direct_gaussian_measurement_covariance_estimator.h"
#include "mixed_generic_direct_gaussian_measurement_covariance_estimator_output.h"
#include "utils.h"

MixedGenericDirectGaussianMeasurementCovarianceEstimator::MixedGenericDirectGaussianMeasurementCovarianceEstimator()
  :MeasurementCovarianceEstimator()
{
}

MixedGenericDirectGaussianMeasurementCovarianceEstimator::MixedGenericDirectGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                                             size_t* seed_in,
                                                                             Data* data_in,
                                                                                                                   std::shared_ptr<Transform> transform_in,
                                                                                                                   std::shared_ptr<Transform> summary_statistics_in,
                                                                                                                   Data* prior_data_in,
                                                                             SimulateModelPtr simulator_in,
                                                                                                                   const std::vector<std::string> &prior_measurement_variables_in,
                                                                                                                   const std::vector<arma::mat> &prior_measurement_noises_in)
: MeasurementCovarianceEstimator(rng_in,
                                 seed_in,
                                 data_in,
                                 transform_in,
                                 summary_statistics_in)
{
  this->set_using_parameters = false;
  this->prior_measurement_variables = prior_measurement_variables_in;
  for (size_t i=0;
       i<this->prior_measurement_variables.size();
       ++i)
  {
    this->kernel.set_mean(this->prior_measurement_variables[i],
                          arma::colvec(prior_measurement_noises_in[i].n_rows));
    this->kernel.set_covariance(this->prior_measurement_variables[i],
                                prior_measurement_noises_in[i]);
  }
  this->transform_function = NULL;
  
  this->simulator = simulator_in;
  this->prior_data = prior_data_in;
}

MixedGenericDirectGaussianMeasurementCovarianceEstimator::~MixedGenericDirectGaussianMeasurementCovarianceEstimator()
{
}

MixedGenericDirectGaussianMeasurementCovarianceEstimator::MixedGenericDirectGaussianMeasurementCovarianceEstimator(const MixedGenericDirectGaussianMeasurementCovarianceEstimator &another)
  :MeasurementCovarianceEstimator(another)
{
  this->make_copy(another);
}

void MixedGenericDirectGaussianMeasurementCovarianceEstimator::operator=(const MixedGenericDirectGaussianMeasurementCovarianceEstimator &another)
{
  if(this == &another)
    return;

  MeasurementCovarianceEstimator::operator=(another);
  this->make_copy(another);
}

void MixedGenericDirectGaussianMeasurementCovarianceEstimator::make_copy(const MixedGenericDirectGaussianMeasurementCovarianceEstimator &another)
{
  //this->conditioned_on_parameters = another.conditioned_on_parameters;
  //this->measurement_kernel = another.measurement_kernel;
  //this->measurement_noise = another.measurement_noise;
  //this->measurement_noise_function = another.measurement_noise_function;
  //this->measurement_kernel_function = another.measurement_kernel_function;
  this->measurement_dimension = another.measurement_dimension;
  this->simulator = another.simulator;
  this->Cygivenx = another.Cygivenx;
  this->likelihood_Cygivenx = another.likelihood_Cygivenx;
  this->gaussian_simulator = another.gaussian_simulator;
  
  this->kernel = another.kernel;
  this->prior_measurement_noise_functions = another.prior_measurement_noise_functions;
  this->transform_function = another.transform_function;
  
  this->prior_measurement_variables = another.prior_measurement_variables;
  this->prior_data = another.prior_data;
}

MeasurementCovarianceEstimator* MixedGenericDirectGaussianMeasurementCovarianceEstimator::duplicate() const
{
  return( new MixedGenericDirectGaussianMeasurementCovarianceEstimator(*this));
}

MeasurementCovarianceEstimatorOutput* MixedGenericDirectGaussianMeasurementCovarianceEstimator::initialise_measurement_covariance_estimator()
{
  MeasurementCovarianceEstimatorOutput* output = new MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput(this);
  return output;
}

MeasurementCovarianceEstimatorOutput* MixedGenericDirectGaussianMeasurementCovarianceEstimator::initialise_measurement_covariance_estimator(const Parameters &conditioned_on_parameters)
{
  MeasurementCovarianceEstimatorOutput* output = new MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput(this);
  return output;
}

void MixedGenericDirectGaussianMeasurementCovarianceEstimator::setup()
{
  this->setup_measurement_variables();
}

void MixedGenericDirectGaussianMeasurementCovarianceEstimator::setup(const Parameters &parameters)
{
  this->setup_measurement_variables(parameters);
}

void MixedGenericDirectGaussianMeasurementCovarianceEstimator::setup_measurement_variables()
{
  Data likelihood_dummy_data = this->simulator(*this->rng,Parameters());
  this->measurement_variables = likelihood_dummy_data.get_vector_variables();
  std::vector<arma::colvec> means;
  means.reserve(this->measurement_variables.size());
  std::vector<arma::mat> covs;
  covs.reserve(this->measurement_variables.size());
  this->measurement_dimension = 0;
  for (size_t i=0; i<this->measurement_variables.size(); ++i)
  {
    arma::colvec current_vector = arma::vectorise(likelihood_dummy_data[this->measurement_variables[i]]);
    means.push_back(arma::zeros<arma::colvec>(current_vector.n_rows));
    covs.push_back(arma::eye(current_vector.n_rows,current_vector.n_rows));
    this->measurement_dimension = this->measurement_dimension + current_vector.n_rows;
  }
  this->gaussian_simulator = GaussianIndependentProposalKernel(this->measurement_variables,
                                                               means,
                                                               covs);
}

void MixedGenericDirectGaussianMeasurementCovarianceEstimator::setup_measurement_variables(const Parameters &conditioned_on_parameters)
{
  Data likelihood_dummy_data = this->simulator(*this->rng,conditioned_on_parameters);
  this->measurement_variables = likelihood_dummy_data.get_vector_variables();
  std::vector<arma::colvec> means;
  means.reserve(this->measurement_variables.size());
  std::vector<arma::mat> covs;
  covs.reserve(this->measurement_variables.size());
  this->measurement_dimension = 0;
  for (size_t i=0; i<this->measurement_variables.size(); ++i)
  {
    arma::colvec current_vector = arma::vectorise(likelihood_dummy_data[this->measurement_variables[i]]);
    means.push_back(arma::zeros<arma::colvec>(current_vector.n_rows));
    covs.push_back(arma::eye(current_vector.n_rows,current_vector.n_rows));
    this->measurement_dimension = this->measurement_dimension + current_vector.n_rows;
  }
  this->gaussian_simulator = GaussianIndependentProposalKernel(this->measurement_variables,
                                                               means,
                                                               covs);
}

/*
arma::mat MixedGenericDirectGaussianMeasurementCovarianceEstimator::get_measurement_covariance() const
{
  arma::mat measurement_noise;
  measurement_noise.zeros(this->measurement_dimension,this->measurement_dimension);
  return measurement_noise;
}
*/

arma::mat MixedGenericDirectGaussianMeasurementCovarianceEstimator::get_Cygivenx() const
{
  return this->Cygivenx;
}

arma::mat MixedGenericDirectGaussianMeasurementCovarianceEstimator::get_adjustment(const arma::mat &Zf,
                                                                                   const arma::mat &Dhathalf,
                                                                                   const arma::mat &P,
                                                                                   const arma::mat &Vtranspose,
                                                                                   const arma::mat &Yhat,
                                                                                   double inverse_incremental_temperature)
{
  // follows https://arxiv.org/abs/2006.02941
  //arma::mat for_eig = V*((inverse_incremental_temperature-1.0)*this->Cygivenx + inverse_incremental_temperature*this->get_prior_measurement_covariance_embedded_in_full_space())*V.t();
  
  arma::mat I;
  I.eye(Vtranspose.n_cols,Vtranspose.n_cols);
  
  arma::mat for_eig = Vtranspose*(arma::inv_sympd(I + Yhat*arma::inv_sympd((inverse_incremental_temperature-1.0)*this->get_Cygivenx() + inverse_incremental_temperature*this->get_prior_measurement_covariance_embedded_in_full_space())*Yhat.t()))*Vtranspose.t();
  
  //std::cout << for_eig << std::endl;
  
  arma::mat U;
  arma::vec diagD;
  //arma::mat Utrans;
  arma::eig_sym(diagD,U,for_eig);
  
  //std::cout << U << std::endl;
  
  arma::mat Dsqrt(diagD.n_elem,diagD.n_elem);
  Dsqrt.diag() = arma::sqrt(diagD);

  /*
  std::cout << P.t() << std::endl;
  std::cout << arma::pinv(Dhathalf) << std::endl;
  std::cout << arma::pinv(Dhathalf)*P.t() << std::endl;
  std::cout << Dsqrt << std::endl;
  std::cout << Dsqrt*arma::pinv(Dhathalf)*P.t() << std::endl;
  std::cout << U << std::endl;
  std::cout << U*Dsqrt*arma::pinv(Dhathalf)*P.t() << std::endl;
  std::cout << Dhathalf << std::endl;
  std::cout << Dhathalf*U*Dsqrt*arma::pinv(Dhathalf)*P.t() << std::endl;
  std::cout << P << std::endl;
  std::cout << P*Dhathalf*U*Dsqrt*arma::pinv(Dhathalf)*P.t() << std::endl;
  */
  
  return P*Dhathalf*U*Dsqrt*arma::pinv(Dhathalf)*P.t();
}

/*
arma::mat MixedGenericDirectGaussianMeasurementCovarianceEstimator::get_adjustment(const arma::mat &Zf,
                                                                                   const arma::mat &Ginv,
                                                                                   const arma::mat &Ftranspose,
                                                                                   const arma::mat &V,
                                                                                   double inverse_incremental_temperature)
{
  // follows https://arxiv.org/abs/2006.02941
  arma::mat for_eig = V*((inverse_incremental_temperature-1.0)*this->Cygivenx + inverse_incremental_temperature*this->get_prior_measurement_covariance_embedded_in_full_space())*V.t();
  
  arma::mat C;
  arma::vec diagGamma;
  arma::mat Ctrans;
  arma::svd(C,diagGamma,Ctrans,for_eig);
  
  arma::mat Gamma(diagGamma.n_elem,diagGamma.n_elem);
  Gamma.diag() = diagGamma;
  arma::mat I;
  I.eye(diagGamma.n_elem,diagGamma.n_elem);
  
  return Zf*C*arma::sqrtmat_sympd(arma::inv_sympd(I+Gamma))*Ginv*Ftranspose;
}
*/

void MixedGenericDirectGaussianMeasurementCovarianceEstimator::set_parameters(const Parameters &conditioned_on_parameters_in)
{
  if (this->set_using_parameters)
  {
    //this->conditioned_on_parameters = conditioned_on_parameters_in;
    for (size_t i=0;
         i<this->prior_measurement_variables.size();
         ++i)
    {
      this->kernel.set_covariance(this->prior_measurement_variables[i],
                                  this->prior_measurement_noise_functions[i](conditioned_on_parameters_in));
    }
  }
}

bool MixedGenericDirectGaussianMeasurementCovarianceEstimator::need_Cxx() const
{
  return true;
}

/*
void MixedGenericDirectGaussianMeasurementCovarianceEstimator::find_partial_Cygivenx(const arma::mat &Cxy,
                                                                                     const arma::mat &Cyy,
                                                                                     const arma::mat &packed_members)
{
  
}
*/

void MixedGenericDirectGaussianMeasurementCovarianceEstimator::find_Cygivenx(const arma::mat &inv_Cxx,
                                                                             const arma::mat &Cxy,
                                                                             const arma::mat &Cyy)
{
  arma::mat subCxy = Cxy.submat(0, 0, Cxy.n_rows-1, this->measurement_dimension-1);
  this->Cygivenx = arma::mat(Cyy.n_rows,Cyy.n_cols);
  this->likelihood_Cygivenx = Cyy.submat(0, 0, this->measurement_dimension-1, this->measurement_dimension-1) - subCxy.t()*inv_Cxx*subCxy;
  this->Cygivenx.submat(0, 0, this->measurement_dimension-1, this->measurement_dimension-1) = this->likelihood_Cygivenx;
}

arma::mat MixedGenericDirectGaussianMeasurementCovarianceEstimator::get_unconditional_measurement_covariance(const arma::mat &Cyy,
                                                                                                             double inverse_incremental_temperature) const
{
  arma::mat gaussian_meas_cov = this->get_prior_measurement_covariance();
  arma::mat total_gaussian_meas_cov = arma::mat(Cyy.n_rows,Cyy.n_cols);
  total_gaussian_meas_cov.submat(this->measurement_dimension, this->measurement_dimension, this->measurement_dimension+gaussian_meas_cov.n_rows-1, this->measurement_dimension+gaussian_meas_cov.n_rows-1) = gaussian_meas_cov;
  return Cyy + (inverse_incremental_temperature-1.0)*this->Cygivenx + inverse_incremental_temperature*total_gaussian_meas_cov;
}

arma::mat MixedGenericDirectGaussianMeasurementCovarianceEstimator::get_prior_measurement_covariance() const
{
  return this->kernel.get_covariance(this->prior_measurement_variables);
}

arma::mat MixedGenericDirectGaussianMeasurementCovarianceEstimator::get_prior_measurement_covariance_embedded_in_full_space() const
{
  arma::mat output(this->Cygivenx.n_rows,this->Cygivenx.n_rows);
  output.submat(this->measurement_dimension, this->measurement_dimension, this->Cygivenx.n_rows-1, this->Cygivenx.n_rows-1) = this->get_prior_measurement_covariance();
  return output;
}

arma::mat MixedGenericDirectGaussianMeasurementCovarianceEstimator::get_measurement_covariance_for_likelihood_ratio(double inverse_incremental_temperature) const
{
  arma::mat likelihood_part = inverse_incremental_temperature*this->Cygivenx;
  arma::mat prior_part = inverse_incremental_temperature*this->get_prior_measurement_covariance();
  arma::mat output(likelihood_part.n_rows+prior_part.n_rows,likelihood_part.n_rows+prior_part.n_rows);
  output.submat(0, 0, likelihood_part.n_rows-1, likelihood_part.n_rows-1) = likelihood_part;
  output.submat(likelihood_part.n_rows, likelihood_part.n_rows, likelihood_part.n_rows+prior_part.n_rows-1, likelihood_part.n_rows+prior_part.n_rows-1) = prior_part;
  return output;
}

Parameters MixedGenericDirectGaussianMeasurementCovarianceEstimator::simulate(const Parameters &current_state)
{
  return this->simulator(*this->rng,
                         current_state);
}

arma::colvec MixedGenericDirectGaussianMeasurementCovarianceEstimator::gaussian_simulate()
{
  return this->gaussian_simulator.independent_simulate(*this->rng).get_colvec(this->measurement_variables);
}

void MixedGenericDirectGaussianMeasurementCovarianceEstimator::change_data()
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
  else
  {
    Rcpp::stop("MixedGenericDirectGaussianMeasurementCovarianceEstimator::change_data - need measurement_variables.");
  }
  
  if (this->prior_measurement_variables.size()>0)
  {
    this->measurement = join_cols(this->measurement,(*this->prior_data)[this->prior_measurement_variables[0]].as_col());
    
    for (size_t i=1; i<this->prior_measurement_variables.size(); ++i)
    {
      this->measurement = join_cols(this->measurement,(*this->prior_data)[this->prior_measurement_variables[i]].as_col());
    }
  }
  else
  {
    Rcpp::stop("MixedGenericDirectGaussianMeasurementCovarianceEstimator::change_data - need prior_measurement_variables.");
  }
}

void MixedGenericDirectGaussianMeasurementCovarianceEstimator::change_data(Data* new_data)
{
  this->current_data = new_data;
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
  
  if (this->prior_measurement_variables.size()>0)
  {
    this->measurement = (*this->prior_data)[this->prior_measurement_variables[0]].as_col();
    
    for (size_t i=1; i<this->prior_measurement_variables.size(); ++i)
    {
      this->measurement = join_cols(this->measurement,(*this->prior_data)[this->prior_measurement_variables[i]].as_col());
    }
  }
}

void MixedGenericDirectGaussianMeasurementCovarianceEstimator::precompute_gaussian_covariance(double inverse_incremental_temperature)
{
  arma::mat for_precomp = this->get_measurement_covariance_for_likelihood_ratio(inverse_incremental_temperature);
  this->inv_sigma_precomp = arma::inv_sympd(for_precomp);
  this->log_det_precomp = arma::log_det_sympd(for_precomp);
}
