#include "density_likelihood_estimator_output.h"
#include "density_likelihood_estimator.h"
#include "density_estimator.h"
#include "density_likelihood_estimator_worker.h"
#include "filesystem.h"
#include "independent_proposal_kernel.h"
#include "density_estimator_output.h"

DensityLikelihoodEstimatorOutput::DensityLikelihoodEstimatorOutput()
  :LikelihoodEstimatorOutput()
{
  this->log_likelihood_smcfixed_part = 0.0;
  this->subsample_log_likelihood_smcfixed_part = 0.0;
}

DensityLikelihoodEstimatorOutput::DensityLikelihoodEstimatorOutput(DensityLikelihoodEstimator* estimator_in,
                                                                   
                                                                   DensityEstimator* density_estimator_in,
                                                                   DensityEstimator* subsample_density_estimator_in)
  :LikelihoodEstimatorOutput()
{
  this->estimator = estimator_in;
  
  if (density_estimator_in!=NULL)
    this->density_estimator_output = density_estimator_in->initialise();
  else
    this->density_estimator_output = NULL;
  
  if (subsample_density_estimator_in!=NULL)
    this->subsample_density_estimator_output = subsample_density_estimator_in->initialise();
  else
    this->subsample_density_estimator_output = NULL;
  //this->density_estimator = this->estimator->density_estimator->duplicate();
  
  this->log_likelihood_smcfixed_part = 0.0;
  this->subsample_log_likelihood_smcfixed_part = 0.0;
  //this->results_name = results_name_in;
}

/*
DensityLikelihoodEstimatorOutput::DensityLikelihoodEstimatorOutput(DensityLikelihoodEstimator* estimator_in,
                                                                   const Parameters &conditioned_on_parameters)
:LikelihoodEstimatorOutput()
{
  this->estimator = estimator_in;
  //this->density_estimator = this->estimator->density_estimator->duplicate();
  this->log_likelihood_smcfixed_part = 0.0;
  this->subsample_log_likelihood_smcfixed_part = 0.0;
  //this->results_name = results_name_in;
}
*/

DensityLikelihoodEstimatorOutput::~DensityLikelihoodEstimatorOutput()
{
  //if (this->density_estimator!=NULL)
  //  delete this->density_estimator;
}

//Copy constructor for the DensityLikelihoodEstimatorOutput class.
DensityLikelihoodEstimatorOutput::DensityLikelihoodEstimatorOutput(const DensityLikelihoodEstimatorOutput &another)
  :LikelihoodEstimatorOutput(another)
{
  this->make_copy(another);
}

void DensityLikelihoodEstimatorOutput::operator=(const DensityLikelihoodEstimatorOutput &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  //if (this->density_estimator!=NULL)
  //  delete this->density_estimator;

  LikelihoodEstimatorOutput::operator=(another);
  this->make_copy(another);
}

LikelihoodEstimatorOutput* DensityLikelihoodEstimatorOutput::duplicate() const
{
  return( new DensityLikelihoodEstimatorOutput(*this));
}

void DensityLikelihoodEstimatorOutput::make_copy(const DensityLikelihoodEstimatorOutput &another)
{
  if (another.density_estimator_output!=NULL)
    this->density_estimator_output = another.density_estimator_output->duplicate();
  else
    this->density_estimator_output = NULL;
  
  if (another.subsample_density_estimator_output!=NULL)
    this->subsample_density_estimator_output = another.subsample_density_estimator_output->duplicate();
  else
    this->subsample_density_estimator_output = NULL;
  
  this->estimator = another.estimator;
  this->log_likelihood_smcfixed_part = another.log_likelihood_smcfixed_part;
  this->subsample_log_likelihood_smcfixed_part = another.subsample_log_likelihood_smcfixed_part;
  //if (another.density_estimator!=NULL)
  //  this->density_estimator = another.density_estimator->duplicate();
  this->points = another.points;
  //this->results_name = another.results_name;
  //this->variables = another.variables;
}

void DensityLikelihoodEstimatorOutput::simulate()
{
  this->estimator->the_worker->simulate(this);
  //Particles particles(this->the_worker->get_particles());
  //this->points = this->estimator->the_worker->get_points();
  //this->density_estimator_output->fit(this->points);
}

void DensityLikelihoodEstimatorOutput::simulate(const Parameters &conditioned_on_parameters)
{
  this->estimator->the_worker->simulate(this,
                                        conditioned_on_parameters);
  //Particles particles(this->the_worker->get_particles());
  //this->points = this->estimator->the_worker->get_points();
  //this->density_estimator_output->fit(this->points);
}

void DensityLikelihoodEstimatorOutput::evaluate_smcfixed_part(const Parameters &parameters)
{
  if (this->estimator->smcfixed_flag)
    this->log_likelihood_smcfixed_part = this->density_estimator_output->evaluate(*this->estimator->data);
}

void DensityLikelihoodEstimatorOutput::evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters)
{
  if (this->estimator->smcfixed_flag)
    this->log_likelihood = this->log_likelihood_smcfixed_part;
  else
    this->log_likelihood = this->density_estimator_output->evaluate(*this->estimator->data);
}

void DensityLikelihoodEstimatorOutput::subsample_simulate(const Parameters &conditioned_on_parameters)
{
  this->estimator->the_worker->subsample_simulate(this,
                                                  conditioned_on_parameters);
  //Particles particles(this->the_worker->get_particles());
  //this->points = this->estimator->the_worker->get_points();
  this->subsample_density_estimator_output->fit(this->points);
}

void DensityLikelihoodEstimatorOutput::subsample_evaluate_smcfixed_part(const Parameters &parameters)
{
  if (this->estimator->smcfixed_flag)
    this->subsample_log_likelihood_smcfixed_part = this->subsample_density_estimator_output->evaluate(*this->estimator->data);
}

void DensityLikelihoodEstimatorOutput::subsample_evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters)
{
  if (this->estimator->smcfixed_flag)
    this->subsample_log_likelihood = this->subsample_log_likelihood_smcfixed_part;
  else
    this->subsample_log_likelihood = this->subsample_density_estimator_output->evaluate(*this->estimator->data);
}

LikelihoodEstimator* DensityLikelihoodEstimatorOutput::get_likelihood_estimator() const
{
  return this->estimator;
}

arma::mat DensityLikelihoodEstimatorOutput::get_gradient_of_log(const std::string &variable,
                                                                const Parameters &x)
{
  // need to use chain rule
  // need differential of llhd wrt parameters (easy) and differential of parameters wrt to theta (not easy)
  Rcpp::stop("DensityLikelihoodEstimatorOutput::get_gradient_of_log - not yet implemented.");
}

arma::mat DensityLikelihoodEstimatorOutput::subsample_get_gradient_of_log(const std::string &variable,
                                                                const Parameters &x)
{
  // need to use chain rule
  // need differential of llhd wrt parameters (easy) and differential of parameters wrt to theta (not easy)
  Rcpp::stop("DensityLikelihoodEstimatorOutput::subsample_get_gradient_of_log - not yet implemented.");
}

void DensityLikelihoodEstimatorOutput::fit(const std::vector<Parameters> &points)
{
  if (this->density_estimator_output!=NULL)
    this->density_estimator_output->fit(points);
}

void DensityLikelihoodEstimatorOutput::subsample_fit(const std::vector<Parameters> &points)
{
  if (this->subsample_density_estimator_output!=NULL)
    this->subsample_density_estimator_output->fit(points);
}

void DensityLikelihoodEstimatorOutput::print(std::ostream &os) const
{

}

void DensityLikelihoodEstimatorOutput::write_to_file(const std::string &dir_name,
                                                     const std::string &index)
{
  std::string directory_name = dir_name + "_density";
  
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
  
  if (!this->estimator->file_stream.is_open())
  {
    this->estimator->file_stream.open(directory_name + "/vector_points.txt");
  }
  if (this->estimator->file_stream.is_open())
  {
    for (auto i = this->points.begin();
         i!=this->points.end();
         ++i)
    {
      this->estimator->file_stream << i->get_rowvec(this->estimator->variables);
    }
    
    //file_stream.close();
  }
  else
  {
    Rcpp::stop("File " + directory_name + "/vector_points.txt" + "cannot be opened.");
  }
}

void DensityLikelihoodEstimatorOutput::forget_you_were_already_written_to_file()
{
}

void DensityLikelihoodEstimatorOutput::close_ofstreams()
{
  this->estimator->log_likelihood_file_stream.close();
  this->estimator->file_stream.close();
}
