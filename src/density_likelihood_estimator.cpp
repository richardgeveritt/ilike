#include "density_likelihood_estimator.h"
#include "density_likelihood_estimator_output.h"
#include "density_likelihood_estimator_worker.h"
#include "density_estimator.h"
#include "parameter_particle_simulator.h"
#include "sequential_density_likelihood_estimator_worker.h"

DensityLikelihoodEstimator::DensityLikelihoodEstimator()
  :LikelihoodEstimator()
{
}

DensityLikelihoodEstimator::DensityLikelihoodEstimator(RandomNumberGenerator* rng_in,
                                                       size_t* seed_in,
                                                       Data* data_in,
                                                       size_t number_of_points_in,
                                                       bool smcfixed_flag_in,
                                                       DensityEstimator* density_estimator_in,
                                                       SimulateIndependentProposalPtr simulate_distribution_in,
                                                       bool parallel_in)
:LikelihoodEstimator(rng_in, seed_in, data_in)
{
  this->density_estimator = density_estimator_in;
  this->subsample_density_estimator = density_estimator_in->duplicate();
  this->simulate_distribution = simulate_distribution_in;
  this->smcfixed_flag = smcfixed_flag_in;
  this->number_of_points = number_of_points_in;

  if (parallel_in==TRUE)
  {
      //this->the_worker = new RcppParallelDensityLikelihoodEstimatorWorker(this);
  }
  else
  {
    this->the_worker = new SequentialDensityLikelihoodEstimatorWorker(this);
  }
}

DensityLikelihoodEstimator::~DensityLikelihoodEstimator()
{
  if (this->density_estimator!=NULL)
    delete this->density_estimator;
  
  if (this->subsample_density_estimator!=NULL)
    delete this->subsample_density_estimator;
  
  if (this->the_worker!=NULL)
    delete this->the_worker;
}

//Copy constructor for the DensityLikelihoodEstimator class.
DensityLikelihoodEstimator::DensityLikelihoodEstimator(const DensityLikelihoodEstimator &another)
  :LikelihoodEstimator(another)
{
  this->make_copy(another);
}

void DensityLikelihoodEstimator::operator=(const DensityLikelihoodEstimator &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  if (this->density_estimator!=NULL)
    delete this->density_estimator;
  
  if (this->subsample_density_estimator!=NULL)
    delete this->subsample_density_estimator;
  
  if (this->the_worker!=NULL)
    delete this->the_worker;
  
  LikelihoodEstimator::operator=(another);
  this->make_copy(another);
}

LikelihoodEstimator* DensityLikelihoodEstimator::duplicate(void)const
{
  return( new DensityLikelihoodEstimator(*this));
}

void DensityLikelihoodEstimator::make_copy(const DensityLikelihoodEstimator &another)
{
  if (another.density_estimator!=NULL)
    this->density_estimator = another.density_estimator->duplicate();
  
  if (another.subsample_density_estimator!=NULL)
    this->subsample_density_estimator = another.subsample_density_estimator->duplicate();
  
  if (another.the_worker!=NULL)
    this->the_worker = another.the_worker->duplicate();
  
  this->smcfixed_flag = another.smcfixed_flag;
  
  this->simulate_distribution = another.simulate_distribution;
  
  this->subsample_simulate_distribution = another.subsample_simulate_distribution;
  
  this->number_of_points = another.number_of_points;
}

// double DensityLikelihoodEstimator::estimate_log_likelihood(const List &inputs,
//                                                          const List &auxiliary_variables) const
// {
//   return this->func(inputs,this->observed_data);
// }

LikelihoodEstimatorOutput* DensityLikelihoodEstimator::initialise()
{
  return new DensityLikelihoodEstimatorOutput(this);
}

LikelihoodEstimatorOutput* DensityLikelihoodEstimator::initialise(const Parameters &parameters)
{
  return new DensityLikelihoodEstimatorOutput(this);
}

//double DensityLikelihoodEstimator::evaluate(const Parameters &parameters)
//{
//  return this->func(parameters,*this->data);
//}

// void DensityLikelihoodEstimator::is_setup_likelihood_estimator(const std::vector<List> &all_points,
//                                                              const std::vector<List> &all_auxiliary_variables)
// {
//
// }
