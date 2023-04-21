#include "exact_likelihood_estimator_output.h"
#include "exact_likelihood_estimator.h"
#include "utils.h"
#include "filesystem.h"

ExactLikelihoodEstimatorOutput::ExactLikelihoodEstimatorOutput()
  :LikelihoodEstimatorOutput()
{
  this->log_likelihood_smcfixed_part = 0.0;
  this->subsample_log_likelihood_smcfixed_part = 0.0;
}

ExactLikelihoodEstimatorOutput::ExactLikelihoodEstimatorOutput(ExactLikelihoodEstimator* estimator_in)
  :LikelihoodEstimatorOutput()
{
  this->estimator = estimator_in;
  this->log_likelihood_smcfixed_part = 0.0;
  this->subsample_log_likelihood_smcfixed_part = 0.0;
}

ExactLikelihoodEstimatorOutput::~ExactLikelihoodEstimatorOutput()
{

}

//Copy constructor for the ExactLikelihoodEstimatorOutput class.
ExactLikelihoodEstimatorOutput::ExactLikelihoodEstimatorOutput(const ExactLikelihoodEstimatorOutput &another)
  :LikelihoodEstimatorOutput(another)
{
  this->make_copy(another);
}

void ExactLikelihoodEstimatorOutput::operator=(const ExactLikelihoodEstimatorOutput &another)
{
  if(this == &another){ //if a==a
    return;
  }

  LikelihoodEstimatorOutput::operator=(another);
  this->make_copy(another);
}

LikelihoodEstimatorOutput* ExactLikelihoodEstimatorOutput::duplicate() const
{
  return( new ExactLikelihoodEstimatorOutput(*this));
}

void ExactLikelihoodEstimatorOutput::make_copy(const ExactLikelihoodEstimatorOutput &another)
{
  this->estimator = another.estimator;
  this->log_likelihood_smcfixed_part = another.log_likelihood_smcfixed_part;
  this->subsample_log_likelihood_smcfixed_part = another.subsample_log_likelihood_smcfixed_part;
}

void ExactLikelihoodEstimatorOutput::simulate()
{
}

void ExactLikelihoodEstimatorOutput::simulate(const Parameters &parameters)
{
}

void ExactLikelihoodEstimatorOutput::evaluate_smcfixed_part(const Parameters &parameters)
{
  if (this->estimator->smcfixed_flag)
    this->log_likelihood_smcfixed_part = this->estimator->evaluate(parameters);
}

void ExactLikelihoodEstimatorOutput::evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters)
{
  if (this->estimator->smcfixed_flag)
    this->log_likelihood = this->log_likelihood_smcfixed_part;
  else
    this->log_likelihood = this->estimator->evaluate(parameters);
}

void ExactLikelihoodEstimatorOutput::subsample_simulate(const Parameters &parameters)
{
}

void ExactLikelihoodEstimatorOutput::subsample_evaluate_smcfixed_part(const Parameters &parameters)
{
  if (this->estimator->smcfixed_flag)
    this->subsample_log_likelihood_smcfixed_part = this->estimator->subsample_evaluate(parameters);
}

void ExactLikelihoodEstimatorOutput::subsample_evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters)
{
  if (this->estimator->smcfixed_flag)
    this->subsample_log_likelihood = this->subsample_log_likelihood_smcfixed_part;
  else
    this->subsample_log_likelihood = this->estimator->subsample_evaluate(parameters);
}

LikelihoodEstimator* ExactLikelihoodEstimatorOutput::get_likelihood_estimator() const
{
  return this->estimator;
}

arma::mat ExactLikelihoodEstimatorOutput::get_gradient_of_log(const std::string &variable,
                                                              const Parameters &x)
{
  return this->estimator->evaluate_gradient(variable,x);
}

arma::mat ExactLikelihoodEstimatorOutput::subsample_get_gradient_of_log(const std::string &variable,
                                                                        const Parameters &x)
{
  return this->estimator->subsample_evaluate_gradient(variable,x);
}

void ExactLikelihoodEstimatorOutput::print(std::ostream &os) const
{

}

void ExactLikelihoodEstimatorOutput::write_to_file(const std::string &dir_name,
                                                   const std::string &index)
{
  std::string directory_name = dir_name + "_exact";
  
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
    this->estimator->log_likelihood_file_stream << this->log_likelihood << std::endl;
    //log_likelihood_file_stream.close();
  }
  else
  {
    Rcpp::stop("File " + directory_name + "/log_likelihood.txt" + "cannot be opened.");
  }
}

void ExactLikelihoodEstimatorOutput::close_ofstreams()
{
  this->estimator->log_likelihood_file_stream.close();
}
