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
  this->estimator_output.reset(estimator_output_in);
}

AnnealedLikelihoodEstimatorOutput::~AnnealedLikelihoodEstimatorOutput()
{
  this->estimator_output.reset();
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
  
  this->estimator_output.reset();
  
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
    this->estimator_output.reset(another.estimator_output->duplicate());
  else
    this->estimator_output.reset();
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

double AnnealedLikelihoodEstimatorOutput::power(const Parameters &parameters) const
{
  if (this->estimator->use_constant)
    return this->estimator->constant_power;
  return this->estimator->function_power(parameters, this->estimator->power_variable);
}

bool AnnealedLikelihoodEstimatorOutput::change_in_place(const Parameters &parameters,
                                                        const std::string &variable,
                                                        size_t row,
                                                        size_t column)
{
  if ((this->estimator_output == NULL) ||
      !this->estimator_output->change_in_place(parameters, variable, row, column))
    return false;
  // (as evaluate_smcadaptive_part_given_smcfixed: a zero power gives zero)
  this->undo_log_likelihood = this->log_likelihood;
  double current_power = this->power(parameters);
  if (current_power != 0.0)
    this->log_likelihood = current_power * this->estimator_output->log_likelihood;
  else
    this->log_likelihood = 0.0;
  return true;
}

void AnnealedLikelihoodEstimatorOutput::undo_in_place()
{
  this->estimator_output->undo_in_place();
  this->log_likelihood = this->undo_log_likelihood;
}

void AnnealedLikelihoodEstimatorOutput::prepare_in_place(const Parameters &parameters)
{
  if (this->estimator_output != NULL)
    this->estimator_output->prepare_in_place(parameters);
}

void AnnealedLikelihoodEstimatorOutput::commit_in_place()
{
  this->estimator_output->commit_in_place();
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

// several variables (C204): the wrapped factor's gradients, times the power
std::vector<arma::mat> AnnealedLikelihoodEstimatorOutput::get_gradients_of_log(const std::vector<std::string> &variables,
                                                                               const Parameters &x)
{
  double power = this->estimator->use_constant ? this->estimator->constant_power
                                               : this->estimator->function_power(x, this->estimator->power_variable);
  std::vector<arma::mat> gradients = this->estimator_output->get_gradients_of_log(variables, x);
  for (arma::mat &gradient : gradients)
  {
    if (power != 0.0)
      gradient = power * gradient;
    else
      gradient = arma::mat(gradient.n_rows, gradient.n_cols, arma::fill::zeros);
  }
  return gradients;
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

// Adaptive summary (dev/plan.md 2.6a): the inner output's value (when the
// inner estimator is smcfixed) or the inner output's own summary; the value
// is then the power at the current target times the inner value, with the
// same rule as evaluate_smcadaptive_part_given_smcfixed for power 0.
int AnnealedLikelihoodEstimatorOutput::adaptive_summary_size() const
{
  if (this->estimator->smcfixed_flag)
    return 1;
  return this->estimator_output->adaptive_summary_size();
}

void AnnealedLikelihoodEstimatorOutput::adaptive_summary(double *values) const
{
  if (this->estimator->smcfixed_flag)
    values[0] = this->estimator_output->log_likelihood;
  else
    this->estimator_output->adaptive_summary(values);
}

double AnnealedLikelihoodEstimatorOutput::smcadaptive_from_summary(const double *values,
                                                                   const Parameters &parameters) const
{
  double inner;
  if (this->estimator->smcfixed_flag)
    inner = values[0];
  else
    inner = this->estimator_output->smcadaptive_from_summary(values, parameters);
  if (this->estimator->use_constant)
  {
    if (this->estimator->constant_power!=0.0)
      return this->estimator->constant_power * inner;
    else
      return 0.0;
  }
  else
  {
    double current_power = this->estimator->function_power(parameters,
                                                           this->estimator->power_variable);
    if (current_power!=0.0)
      return current_power * inner;
    else
      return 0.0;
  }
}

void AnnealedLikelihoodEstimatorOutput::restore_from_summary(const double *values,
                                                             const Parameters &parameters)
{
  // the inner output as its full evaluation left it (summary: its value),
  // then this output's value at the power given by parameters
  this->estimator_output->restore_from_summary(values, parameters);
  this->log_likelihood = this->smcadaptive_from_summary(values, parameters);
}

bool AnnealedLikelihoodEstimatorOutput::repoint_estimator(LikelihoodEstimator* estimator_in)
{
  AnnealedLikelihoodEstimator* annealed = dynamic_cast<AnnealedLikelihoodEstimator*>(estimator_in);
  if (annealed == NULL)
    return false;
  this->estimator = annealed;
  if (this->estimator_output != NULL)
    return this->estimator_output->repoint_estimator(annealed->estimator.get());
  return true;
}

void AnnealedLikelihoodEstimatorOutput::carry_over_from(const LikelihoodEstimatorOutput &previous,
                                                        const Parameters &previous_parameters,
                                                        const Parameters &new_parameters,
                                                        const std::vector<std::string>* changed)
{
  const AnnealedLikelihoodEstimatorOutput* source = dynamic_cast<const AnnealedLikelihoodEstimatorOutput*>(&previous);
  if ((source != NULL) && (source->estimator_output != NULL) && (this->estimator_output != NULL))
    this->estimator_output->carry_over_from(*source->estimator_output, previous_parameters, new_parameters, changed);
}

}
