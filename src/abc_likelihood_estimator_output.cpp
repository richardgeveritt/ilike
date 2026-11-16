#include "abc_likelihood_estimator_output.h"
#include "abc_likelihood_estimator.h"
#include "filesystem.h"
#include "ilike_hdf5_utils.h"

namespace ilike
{

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
  double values[1] = {this->log_likelihood};
  this->write_compact_values(dir_name, index, values);
}

void ABCLikelihoodEstimatorOutput::forget_you_were_already_written_to_file()
{
}

size_t ABCLikelihoodEstimatorOutput::compact_size() const
{
  return 1;
}

void ABCLikelihoodEstimatorOutput::compact_values(double *values) const
{
  values[0] = this->log_likelihood;
}

void ABCLikelihoodEstimatorOutput::write_compact_values(const std::string &dir_name,
                                                        const std::string &index,
                                                        const double *values)
{
  // Buffered: the values of a generation are written in one call by
  // close_ofstreams (this used to stat the directory and extend the HDF5
  // dataset once per point).
  if (this->estimator->h5_file_path.empty())
  {
    std::string directory_name = dir_name + "_abc";
    if (!directory_exists(directory_name))
      make_directory(directory_name);
    this->estimator->h5_file_path = directory_name + "/output.h5";
  }
  this->estimator->pending_log_likelihoods.push_back(values[0]);
}

void ABCLikelihoodEstimatorOutput::close_ofstreams()
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
}


// Adaptive summary (dev/plan.md 2.6a): the distance and scale constant found
// in the smcfixed part; the kernel at the current tolerance is then
// evaluate_kernel_given_distance. When smcfixed, the value itself.
int ABCLikelihoodEstimatorOutput::adaptive_summary_size() const
{
  return this->estimator->smcfixed_flag ? 1 : 2;
}

void ABCLikelihoodEstimatorOutput::adaptive_summary(double *values) const
{
  if (this->estimator->smcfixed_flag)
  {
    values[0] = this->log_likelihood;
  }
  else
  {
    values[0] = this->distance;
    values[1] = this->scale_constant;
  }
}

double ABCLikelihoodEstimatorOutput::smcadaptive_from_summary(const double *values,
                                                              const Parameters &parameters) const
{
  if (this->estimator->smcfixed_flag)
    return values[0];
  return this->estimator->abc_kernel->evaluate_kernel_given_distance(parameters,
                                                                     values[0],
                                                                     values[1]);
}


bool ABCLikelihoodEstimatorOutput::repoint_estimator(LikelihoodEstimator* estimator_in)
{
  ABCLikelihoodEstimator* abc = dynamic_cast<ABCLikelihoodEstimator*>(estimator_in);
  if (abc == NULL)
    return false;
  this->estimator = abc;
  return true;
}


void ABCLikelihoodEstimatorOutput::smcadaptive_from_summaries(const double *values,
                                                              size_t stride,
                                                              size_t n,
                                                              const Parameters &parameters,
                                                              double *results) const
{
  if (this->estimator->smcfixed_flag)
  {
    for (size_t i = 0; i < n; ++i)
      results[i] = values[i*stride];
    return;
  }
  // as smcadaptive_from_summary, with the tolerance looked up once
  double epsilon = this->estimator->abc_kernel->get_epsilon(parameters);
  for (size_t i = 0; i < n; ++i)
    results[i] = this->estimator->abc_kernel->evaluate_kernel_given_distance_at(epsilon,
                                                                                values[i*stride],
                                                                                values[i*stride + 1]);
}

}
