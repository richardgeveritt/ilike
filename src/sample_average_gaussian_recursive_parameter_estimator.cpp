#include "sample_average_gaussian_recursive_parameter_estimator.h"
#include "utils.h"
#include "proposal_kernel.h"

SampleAverageGaussianRecursiveParameterEstimator::SampleAverageGaussianRecursiveParameterEstimator()
  :GaussianRecursiveParameterEstimator()
{
}

SampleAverageGaussianRecursiveParameterEstimator::~SampleAverageGaussianRecursiveParameterEstimator()
{
  
}

//Copy constructor for the SampleAverageGaussianRecursiveParameterEstimator class.
SampleAverageGaussianRecursiveParameterEstimator::SampleAverageGaussianRecursiveParameterEstimator(const SampleAverageGaussianRecursiveParameterEstimator &another)
  :GaussianRecursiveParameterEstimator(another)
{
  this->make_copy(another);
}

void SampleAverageGaussianRecursiveParameterEstimator::operator=(const SampleAverageGaussianRecursiveParameterEstimator &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  GaussianRecursiveParameterEstimator::operator=(another);
  this->make_copy(another);
}

RecursiveParameterEstimator* SampleAverageGaussianRecursiveParameterEstimator::duplicate(void) const
{
  return( new SampleAverageGaussianRecursiveParameterEstimator(*this));
}

GaussianRecursiveParameterEstimator* SampleAverageGaussianRecursiveParameterEstimator::gaussian_duplicate(void) const
{
  return( new SampleAverageGaussianRecursiveParameterEstimator(*this));
}

void SampleAverageGaussianRecursiveParameterEstimator::make_copy(const SampleAverageGaussianRecursiveParameterEstimator &another)
{
  this->gain = another.gain;
}

void SampleAverageGaussianRecursiveParameterEstimator::update(const std::string &variable_name,
                                                              const Particle &latest_particle,
                                                              size_t iteration_counter,
                                                              ProposalKernel* proposal)
{
  if (iteration_counter==0)
  {
    arma::colvec current_point = latest_particle.move_parameters->get_colvec(variable_name);
    this->estimated.get_mean() = current_point;
    this->estimated.get_covariance() = arma::mat(current_point.n_rows,current_point.n_rows,arma::fill::zeros);
  }
  else
  {
    double current_gain = this->gain(iteration_counter+1);
    arma::colvec current_point = latest_particle.move_parameters->get_colvec(variable_name);
    this->estimated.get_mean() = this->estimated.get_mean() + current_gain*(latest_particle.move_parameters->get_colvec(variable_name) - this->estimated.get_mean());
    
    arma::colvec current_difference = current_point - this->estimated.get_mean();
    this->estimated.get_covariance() = this->estimated.get_covariance() + current_gain*(current_difference*current_difference.t() - this->estimated.get_covariance());
  }
}
