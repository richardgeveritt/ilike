#include "abc_likelihood_estimator_output.h"
#include "abc_likelihood_estimator.h"
#include "filesystem.h"

ABCLikelihoodEstimatorOutput::ABCLikelihoodEstimatorOutput()
  :LikelihoodEstimatorOutput()
{
}

ABCLikelihoodEstimatorOutput::ABCLikelihoodEstimatorOutput(ABCLikelihoodEstimator* estimator_in)
  :LikelihoodEstimatorOutput()
{
  this->estimator = estimator_in;
}

ABCLikelihoodEstimatorOutput::~ABCLikelihoodEstimatorOutput()
{
}

//Copy constructor for the ABCLikelihoodEstimatorOutput class.
ABCLikelihoodEstimatorOutput::ABCLikelihoodEstimatorOutput(const ABCLikelihoodEstimatorOutput &another)
  :LikelihoodEstimatorOutput(another)
{
  this->make_copy(another);
}

void ABCLikelihoodEstimatorOutput::operator=(const ABCLikelihoodEstimatorOutput &another)
{
  if(this == &another){ //if a==a
    return;
  }

  LikelihoodEstimatorOutput::operator=(another);
  this->make_copy(another);
}

LikelihoodEstimatorOutput* ABCLikelihoodEstimatorOutput::duplicate() const
{
  return( new ABCLikelihoodEstimatorOutput(*this));
}

void ABCLikelihoodEstimatorOutput::make_copy(const ABCLikelihoodEstimatorOutput &another)
{
  this->estimator = another.estimator;
  this->distance = another.distance;
  this->scale_constant = another.scale_constant;
}

void ABCLikelihoodEstimatorOutput::simulate()
{
}

void ABCLikelihoodEstimatorOutput::simulate(const Parameters &parameters)
{
}

void ABCLikelihoodEstimatorOutput::evaluate_smcfixed_part(const Parameters &parameters)
{
  if (this->estimator->smcfixed_flag)
  {
    this->log_likelihood = this->estimator->abc_kernel->likelihood_evaluate(parameters);
  }
  else
  {
    this->estimator->abc_kernel->find_distance(parameters,
                                               this->distance,
                                               this->scale_constant);
  }
}

void ABCLikelihoodEstimatorOutput::evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters)
{
  if (!this->estimator->smcfixed_flag)
  {
    this->log_likelihood = this->estimator->abc_kernel->evaluate_kernel_given_distance(parameters,
                                                                                       this->distance,
                                                                                       this->scale_constant);
  }
}

void ABCLikelihoodEstimatorOutput::subsample_simulate(const Parameters &parameters)
{
}

void ABCLikelihoodEstimatorOutput::subsample_evaluate_smcfixed_part(const Parameters &parameters)
{
  if (this->estimator->smcfixed_flag)
  {
    this->estimator->abc_kernel->set_data(this->estimator->subsampler->small_data);
    this->subsample_log_likelihood = this->estimator->abc_kernel->likelihood_evaluate(parameters);
  }
  else
  {
    this->estimator->abc_kernel->set_data(this->estimator->subsampler->small_data);
    this->estimator->abc_kernel->find_distance(parameters,
                                               this->distance,
                                               this->scale_constant);
  }
}

void ABCLikelihoodEstimatorOutput::subsample_evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters)
{
  if (!this->estimator->smcfixed_flag)
  {
    this->subsample_log_likelihood = this->estimator->abc_kernel->evaluate_kernel_given_distance(parameters,
                                                                                                 this->distance,
                                                                                                 this->scale_constant);
  }
}

LikelihoodEstimator* ABCLikelihoodEstimatorOutput::get_likelihood_estimator() const
{
  return this->estimator;
}

arma::mat ABCLikelihoodEstimatorOutput::get_gradient_of_log(const std::string &variable,
                                                            const Parameters &x)
{
  Rcpp::stop("ABCLikelihoodEstimatorOutput::get_gradient_of_log - not written yet.");
}

arma::mat ABCLikelihoodEstimatorOutput::subsample_get_gradient_of_log(const std::string &variable,
                                                                 const Parameters &x)
{
  Rcpp::stop("ABCLikelihoodEstimatorOutput::get_gradient_of_log - not written yet.");
}

void ABCLikelihoodEstimatorOutput::print(std::ostream &os) const
{

}

void ABCLikelihoodEstimatorOutput::write_to_file(const std::string &dir_name,
                                                 const std::string &index)
{
  std::string directory_name = dir_name + "_abc";
  
  //if (index!="")
  //  directory_name = directory_name + "_" + index;
  
  if (!directory_exists(directory_name))
  {
    make_directory(directory_name);
  }
  
  if (!this->estimator->log_likelihood_file_stream.is_open())
  {
    this->estimator->log_likelihood_file_stream.open(directory_name + "/log_likelihood.txt",std::ios::out | std::ios::app);
  }
  if (this->estimator->log_likelihood_file_stream.is_open())
  {
    this->estimator->log_likelihood_file_stream << std::setprecision(std::numeric_limits<double>::max_digits10) << this->log_likelihood << std::endl;
    //log_likelihood_file_stream.close();
  }
  else
  {
    Rcpp::stop("File " + directory_name + "/log_likelihood.txt" + "cannot be opened.");
  }
}

void ABCLikelihoodEstimatorOutput::forget_you_were_already_written_to_file()
{
}

void ABCLikelihoodEstimatorOutput::close_ofstreams()
{
  this->estimator->log_likelihood_file_stream.close();
}
