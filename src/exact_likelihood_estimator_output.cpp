#include "exact_likelihood_estimator_output.h"
#include "exact_likelihood_estimator.h"
#include "utils.h"
#include "filesystem.h"
#include "ilike_hdf5_utils.h"

namespace ilike
{
ExactLikelihoodEstimatorOutput::ExactLikelihoodEstimatorOutput()
:LikelihoodEstimatorOutput()
{
  this->estimator = nullptr;
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
  if (!this->estimator) return;
  // Just record the file path on first call; accumulate value in buffer.
  // No HDF5 I/O here — the batch write happens in close_ofstreams().
  if (this->estimator->h5_file_path.empty())
    this->estimator->h5_file_path = dir_name + "_exact.h5";
  
  this->estimator->pending_log_likelihoods.push_back(this->log_likelihood);
}

void ExactLikelihoodEstimatorOutput::forget_you_were_already_written_to_file()
{
}

bool ExactLikelihoodEstimatorOutput::reinitialise(const Parameters &parameters)
{
  // the state of a new ExactLikelihoodEstimatorOutput(estimator) after
  // simulate (nothing to simulate) with write_to_file_flag set to false
  this->log_likelihood = 0.0;
  this->subsample_log_likelihood = 0.0;
  this->write_to_file_flag = false;
  this->log_likelihood_smcfixed_part = 0.0;
  this->subsample_log_likelihood_smcfixed_part = 0.0;
  return true;
}

bool ExactLikelihoodEstimatorOutput::assign_from(const LikelihoodEstimatorOutput &another)
{
  const ExactLikelihoodEstimatorOutput* exact = dynamic_cast<const ExactLikelihoodEstimatorOutput*>(&another);
  if (exact == NULL)
    return false;
  LikelihoodEstimatorOutput::make_copy(*exact);
  this->make_copy(*exact);
  return true;
}

size_t ExactLikelihoodEstimatorOutput::compact_size() const
{
  return 1;
}

void ExactLikelihoodEstimatorOutput::compact_values(double *values) const
{
  values[0] = this->log_likelihood;
}

void ExactLikelihoodEstimatorOutput::write_compact_values(const std::string &dir_name,
                                                          const std::string &index,
                                                          const double *values)
{
  // as write_to_file, with the stored value
  if (!this->estimator) return;
  if (this->estimator->h5_file_path.empty())
    this->estimator->h5_file_path = dir_name + "_exact.h5";
  this->estimator->pending_log_likelihoods.push_back(values[0]);
}

void ExactLikelihoodEstimatorOutput::close_ofstreams()
{
  if (!this->estimator) return;
  // Flush accumulated values in a single HDF5 extend, then clear the buffer.
  if (!this->estimator->pending_log_likelihoods.empty() &&
      !this->estimator->h5_file_path.empty())
  {
    if (!this->estimator->h5_file)
      this->estimator->h5_file = h5_open_or_create(this->estimator->h5_file_path);
    
    h5_append_doubles(*this->estimator->h5_file,
                      "log_likelihood",
                      this->estimator->pending_log_likelihoods);
    this->estimator->pending_log_likelihoods.clear();
  }
  // Close the file and forget its path: the next output written (the next
  // SMC iteration) goes to its own directory. Keeping the path sent every
  // iteration's values to the first iteration's file (B22).
  this->estimator->h5_file.reset();
  this->estimator->h5_file_path.clear();
}

// Adaptive summary (dev/plan.md 2.6a): when the estimator is smcfixed the
// adaptive part is the cached smcfixed value; otherwise the model is
// evaluated again at the new target (no summary).
int ExactLikelihoodEstimatorOutput::adaptive_summary_size() const
{
  return this->estimator->smcfixed_flag ? 1 : -1;
}

void ExactLikelihoodEstimatorOutput::adaptive_summary(double *values) const
{
  values[0] = this->log_likelihood_smcfixed_part;
}

double ExactLikelihoodEstimatorOutput::smcadaptive_from_summary(const double *values,
                                                                const Parameters &parameters) const
{
  return values[0];
}


void ExactLikelihoodEstimatorOutput::restore_from_summary(const double *values,
                                                         const Parameters &parameters)
{
  // after evaluate_smcfixed_part and evaluate_smcadaptive_part_given_smcfixed
  // of an smcfixed estimator: both parts hold the model's value
  this->reinitialise(parameters);
  this->log_likelihood_smcfixed_part = values[0];
  this->log_likelihood = values[0];
}


bool ExactLikelihoodEstimatorOutput::repoint_estimator(LikelihoodEstimator* estimator_in)
{
  ExactLikelihoodEstimator* exact = dynamic_cast<ExactLikelihoodEstimator*>(estimator_in);
  if (exact == NULL)
    return false;
  this->estimator = exact;
  return true;
}

}
