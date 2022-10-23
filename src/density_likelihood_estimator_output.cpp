#include "density_likelihood_estimator_output.h"
#include "density_likelihood_estimator.h"
#include "density_estimator.h"
#include "density_likelihood_estimator_worker.h"

DensityLikelihoodEstimatorOutput::DensityLikelihoodEstimatorOutput()
  :LikelihoodEstimatorOutput()
{
  this->log_likelihood_smcfixed_part = 0.0;
  this->subsample_log_likelihood_smcfixed_part = 0.0;
}

DensityLikelihoodEstimatorOutput::DensityLikelihoodEstimatorOutput(DensityLikelihoodEstimator* estimator_in)
  :LikelihoodEstimatorOutput()
{
  this->estimator = estimator_in;
  //this->density_estimator = this->estimator->density_estimator->duplicate();
  this->log_likelihood_smcfixed_part = 0.0;
  this->subsample_log_likelihood_smcfixed_part = 0.0;
}

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

LikelihoodEstimatorOutput* DensityLikelihoodEstimatorOutput::duplicate(void)const
{
  return( new DensityLikelihoodEstimatorOutput(*this));
}

void DensityLikelihoodEstimatorOutput::make_copy(const DensityLikelihoodEstimatorOutput &another)
{
  this->estimator = another.estimator;
  this->log_likelihood_smcfixed_part = another.log_likelihood_smcfixed_part;
  this->subsample_log_likelihood_smcfixed_part = another.subsample_log_likelihood_smcfixed_part;
  //if (another.density_estimator!=NULL)
  //  this->density_estimator = another.density_estimator->duplicate();
  this->points = another.points;
}

void DensityLikelihoodEstimatorOutput::simulate()
{
  this->estimator->the_worker->simulate();
  //Particles particles(this->the_worker->get_particles());
  this->points = this->estimator->the_worker->get_points();
  this->estimator->density_estimator->fit(this->points);
}

void DensityLikelihoodEstimatorOutput::simulate(const Parameters &conditioned_on_parameters)
{
  this->estimator->the_worker->simulate(conditioned_on_parameters);
  //Particles particles(this->the_worker->get_particles());
  this->points = this->estimator->the_worker->get_points();
  this->estimator->density_estimator->fit(this->points);
}

void DensityLikelihoodEstimatorOutput::evaluate_smcfixed_part(const Parameters &parameters)
{
  if (this->estimator->smcfixed_flag)
    this->log_likelihood_smcfixed_part = this->estimator->density_estimator->evaluate(parameters);
}

void DensityLikelihoodEstimatorOutput::evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters)
{
  if (this->estimator->smcfixed_flag)
    this->log_likelihood = this->log_likelihood_smcfixed_part;
  else
    this->log_likelihood = this->estimator->density_estimator->evaluate(parameters);
}

void DensityLikelihoodEstimatorOutput::subsample_simulate(const Parameters &conditioned_on_parameters)
{
  this->estimator->the_worker->subsample_simulate(conditioned_on_parameters);
  //Particles particles(this->the_worker->get_particles());
  this->points = this->estimator->the_worker->get_points();
  this->estimator->subsample_density_estimator->fit(this->points);
}

void DensityLikelihoodEstimatorOutput::subsample_evaluate_smcfixed_part(const Parameters &parameters)
{
  if (this->estimator->smcfixed_flag)
    this->subsample_log_likelihood_smcfixed_part = this->estimator->subsample_density_estimator->evaluate(parameters);
}

void DensityLikelihoodEstimatorOutput::subsample_evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters)
{
  if (this->estimator->smcfixed_flag)
    this->subsample_log_likelihood = this->subsample_log_likelihood_smcfixed_part;
  else
    this->subsample_log_likelihood = this->estimator->subsample_density_estimator->evaluate(parameters);
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
  throw std::runtime_error("DensityLikelihoodEstimatorOutput::get_gradient_of_log - not yet implemented.");
}

arma::mat DensityLikelihoodEstimatorOutput::subsample_get_gradient_of_log(const std::string &variable,
                                                                const Parameters &x)
{
  // need to use chain rule
  // need differential of llhd wrt parameters (easy) and differential of parameters wrt to theta (not easy)
  throw std::runtime_error("DensityLikelihoodEstimatorOutput::subsample_get_gradient_of_log - not yet implemented.");
}

void DensityLikelihoodEstimatorOutput::print(std::ostream &os) const
{

}
