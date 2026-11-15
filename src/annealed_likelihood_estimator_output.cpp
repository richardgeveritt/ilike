#include "annealed_likelihood_estimator_output.h"
#include "annealed_likelihood_estimator.h"
#include "utils.h"
#include "filesystem.h"
#include "ilike_hdf5_utils.h"

namespace ilike
{
AnnealedLikelihoodEstimatorOutput::AnnealedLikelihoodEstimatorOutput()
:LikelihoodEstimatorOutput()
{
}

AnnealedLikelihoodEstimatorOutput::AnnealedLikelihoodEstimatorOutput(AnnealedLikelihoodEstimator* estimator_in,
                                                                     LikelihoodEstimatorOutput* estimator_output_in)
:LikelihoodEstimatorOutput()
{
  this->estimator = estimator_in;
  this->estimator_output = estimator_output_in;
}

AnnealedLikelihoodEstimatorOutput::~AnnealedLikelihoodEstimatorOutput()
{
  if (this->estimator_output!=NULL)
    delete this->estimator_output;
}

//Copy constructor for the AnnealedLikelihoodEstimatorOutput class.
AnnealedLikelihoodEstimatorOutput::AnnealedLikelihoodEstimatorOutput(const AnnealedLikelihoodEstimatorOutput &another)
:LikelihoodEstimatorOutput(another)
{
  this->make_copy(another);
}

void AnnealedLikelihoodEstimatorOutput::operator=(const AnnealedLikelihoodEstimatorOutput &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  if (this->estimator_output!=NULL)
    delete this->estimator_output;
  
  LikelihoodEstimatorOutput::operator=(another);
  this->make_copy(another);
}

LikelihoodEstimatorOutput* AnnealedLikelihoodEstimatorOutput::duplicate() const
{
  return( new AnnealedLikelihoodEstimatorOutput(*this));
}

void AnnealedLikelihoodEstimatorOutput::make_copy(const AnnealedLikelihoodEstimatorOutput &another)
{
  this->estimator = another.estimator;
  if (another.estimator_output!=NULL)
    this->estimator_output = another.estimator_output->duplicate();
  else
    this->estimator_output = NULL;
}

void AnnealedLikelihoodEstimatorOutput::simulate()
{
  this->estimator_output->simulate();
}

void AnnealedLikelihoodEstimatorOutput::simulate(const Parameters &parameters)
{
  this->estimator_output->simulate(parameters);
}

void AnnealedLikelihoodEstimatorOutput::evaluate_smcfixed_part(const Parameters &parameters)
{
  if (this->estimator->smcfixed_flag)
  {
    this->estimator_output->evaluate(parameters);
  }
  else
  {
    this->estimator_output->evaluate_smcfixed_part(parameters);
  }
}

void AnnealedLikelihoodEstimatorOutput::evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters)
{
  if (!this->estimator->smcfixed_flag)
  {
    this->estimator_output->evaluate_smcadaptive_part_given_smcfixed(parameters);
  }
  
  if (this->estimator->use_constant)
  {
    if (this->estimator->constant_power!=0.0)
      this->log_likelihood = this->estimator->constant_power * this->estimator_output->log_likelihood;
    else
      this->log_likelihood = 0.0;
  }
  else
  {
    double current_power = this->estimator->function_power(parameters,
                                                           this->estimator->power_variable);
    
    if (current_power!=0.0)
      this->log_likelihood = current_power * this->estimator_output->log_likelihood;
    else
      this->log_likelihood = 0.0;
  }
}

void AnnealedLikelihoodEstimatorOutput::subsample_simulate(const Parameters &parameters)
{
  this->estimator_output->subsample_simulate(parameters);
}

void AnnealedLikelihoodEstimatorOutput::subsample_evaluate_smcfixed_part(const Parameters &parameters)
{
  if (this->estimator->smcfixed_flag)
  {
    this->estimator_output->subsample_evaluate(parameters);
  }
  else
  {
    this->estimator_output->subsample_evaluate_smcfixed_part(parameters);
  }
}

void AnnealedLikelihoodEstimatorOutput::subsample_evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters)
{
  if (!this->estimator->smcfixed_flag)
  {
    this->estimator_output->subsample_evaluate_smcadaptive_part_given_smcfixed(parameters);
  }
  
  if (this->estimator->use_constant)
  {
    if (this->estimator->constant_power!=0.0)
    {
      this->subsample_log_likelihood = this->estimator->constant_power * this->estimator_output->subsample_log_likelihood;
    }
    else
    {
      this->subsample_log_likelihood = 0.0;
    }
  }
  else
  {
    double current_power = this->estimator->function_power(parameters,
                                                           this->estimator->power_variable);
    if (current_power!=0.0)
      this->subsample_log_likelihood = current_power * this->estimator_output->subsample_log_likelihood;
    else
      this->subsample_log_likelihood = 0.0;
  }
}

LikelihoodEstimator* AnnealedLikelihoodEstimatorOutput::get_likelihood_estimator() const
{
  return this->estimator;
}

arma::mat AnnealedLikelihoodEstimatorOutput::get_gradient_of_log(const std::string &variable,
                                                                 const Parameters &x)
{
  arma::mat gradient;
  if (this->estimator->use_constant)
  {
    if (this->estimator->constant_power!=0.0)
    {
      gradient = this->estimator->constant_power * this->estimator_output->get_gradient_of_log(variable,
                                                                                               x);
    }
    else
    {
      arma::mat grad = this->estimator_output->get_gradient_of_log(variable,
                                                                   x);
      gradient = arma::mat(grad.n_rows,grad.n_cols);
      gradient.fill(0.0);
    }
  }
  else
  {
    double current_power = this->estimator->function_power(x,
                                                           this->estimator->power_variable);
    if (current_power!=0.0)
    {
      gradient = current_power * this->estimator_output->get_gradient_of_log(variable,
                                                                             x);
    }
    else
    {
      arma::mat grad = this->estimator_output->get_gradient_of_log(variable,
                                                                   x);
      gradient = arma::mat(grad.n_rows,grad.n_cols);
      gradient.fill(0.0);
    }
  }
  return gradient;
}

arma::mat AnnealedLikelihoodEstimatorOutput::subsample_get_gradient_of_log(const std::string &variable,
                                                                           const Parameters &x)
{
  arma::mat gradient;
  if (this->estimator->use_constant)
  {
    if (this->estimator->constant_power!=0.0)
    {
      gradient = this->estimator->constant_power * this->estimator_output->subsample_get_gradient_of_log(variable,
                                                                                                         x);
    }
    else
    {
      arma::mat grad = this->estimator_output->subsample_get_gradient_of_log(variable,
                                                                             x);
      gradient = arma::mat(grad.n_rows,grad.n_cols);
      gradient.fill(0.0);
    }
  }
  else
  {
    double current_power = this->estimator->function_power(x,
                                                           this->estimator->power_variable);
    if (current_power!=0.0)
    {
      gradient = current_power * this->estimator_output->subsample_get_gradient_of_log(variable,
                                                                                       x);
    }
    else
    {
      arma::mat grad = this->estimator_output->subsample_get_gradient_of_log(variable,
                                                                             x);
      gradient = arma::mat(grad.n_rows,grad.n_cols);
      gradient.fill(0.0);
    }
  }
  return gradient;
}

void AnnealedLikelihoodEstimatorOutput::print(std::ostream &os) const
{
  
}

void AnnealedLikelihoodEstimatorOutput::write_to_file(const std::string &dir_name,
                                                      const std::string &index)
{
  // Buffered (see close_ofstreams); the inner output writes as before.
  std::string directory_name = dir_name + "_annealed";
  if (this->estimator->h5_file_path.empty())
  {
    if (!directory_exists(directory_name))
      make_directory(directory_name);
    this->estimator->h5_file_path = directory_name + "/output.h5";
  }
  this->estimator->pending_log_likelihoods.push_back(this->log_likelihood);
  this->estimator_output->write_to_file(directory_name, index);
}

size_t AnnealedLikelihoodEstimatorOutput::compact_size() const
{
  // its own value and the inner output's record, if the inner output has one
  size_t inner = this->estimator_output->compact_size();
  return (inner == 0) ? 0 : 1 + inner;
}

void AnnealedLikelihoodEstimatorOutput::compact_values(double *values) const
{
  values[0] = this->log_likelihood;
  this->estimator_output->compact_values(values + 1);
}

void AnnealedLikelihoodEstimatorOutput::write_compact_values(const std::string &dir_name,
                                                             const std::string &index,
                                                             const double *values)
{
  std::string directory_name = dir_name + "_annealed";
  if (this->estimator->h5_file_path.empty())
  {
    if (!directory_exists(directory_name))
      make_directory(directory_name);
    this->estimator->h5_file_path = directory_name + "/output.h5";
  }
  this->estimator->pending_log_likelihoods.push_back(values[0]);
  this->estimator_output->write_compact_values(directory_name, index, values + 1);
}

void AnnealedLikelihoodEstimatorOutput::forget_you_were_already_written_to_file()
{
  this->estimator_output->forget_you_were_already_written_to_file();
}

void AnnealedLikelihoodEstimatorOutput::close_ofstreams()
{
  if (!this->estimator->pending_log_likelihoods.empty() &&
      !this->estimator->h5_file_path.empty())
  {
    std::shared_ptr<HighFive::File> file = h5_open_or_create(this->estimator->h5_file_path);
    h5_append_doubles(*file, "log_likelihood", this->estimator->pending_log_likelihoods);
    this->estimator->pending_log_likelihoods.clear();
  }
  this->estimator->h5_file.reset();
  this->estimator->h5_file_path.clear();
  this->estimator_output->close_ofstreams();
}
}
