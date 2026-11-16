#include "density_likelihood_estimator.h"
#include "density_likelihood_estimator_output.h"
#include "density_likelihood_estimator_worker.h"
#include "density_estimator.h"
#include "parameter_particle_simulator.h"
#include "sequential_density_likelihood_estimator_worker.h"
#include "independent_proposal_kernel.h"
//#include "transform.h"

namespace ilike
{
DensityLikelihoodEstimator::DensityLikelihoodEstimator()
:LikelihoodEstimator()
{
  //this->transform = NULL;
}

DensityLikelihoodEstimator::DensityLikelihoodEstimator(RandomNumberGenerator* rng_in,
                                                       size_t* seed_in,
                                                       Data* data_in,
                                                       const Parameters &algorithm_parameters_in,
                                                       size_t number_of_points_in,
                                                       bool smcfixed_flag_in,
                                                       DensityEstimator* density_estimator_in,
                                                       IndependentProposalKernel* proposal_in,
                                                       bool make_subsample_version_in,
                                                       bool store_output_in,
                                                       bool parallel_in,
                                                       size_t grain_size_in)
:LikelihoodEstimator(rng_in, seed_in, data_in, algorithm_parameters_in, smcfixed_flag_in)
{
  proposal_in->set_proposal_parameters(&this->algorithm_parameters);
  this->density_estimator.reset(density_estimator_in);
  if (make_subsample_version_in==true)
    this->subsample_density_estimator.reset(density_estimator_in->duplicate());
  else
    this->subsample_density_estimator.reset();
  this->proposal.reset(proposal_in);
  this->subsample_proposal.reset();
  this->number_of_points = number_of_points_in;
  
  this->store_output = store_output_in;
  //this->transform = NULL;
  
  if (parallel_in==TRUE)
  {
    //this->the_worker.reset(new RcppParallelDensityLikelihoodEstimatorWorker(this));
  }
  else
  {
    this->the_worker.reset(new SequentialDensityLikelihoodEstimatorWorker(this));
  }
}

DensityLikelihoodEstimator::~DensityLikelihoodEstimator()
{
  this->density_estimator.reset();
  
  this->subsample_density_estimator.reset();
  
  this->the_worker.reset();
  
  this->proposal.reset();
  
  this->subsample_proposal.reset();
  
  //if (this->transform!=NULL)
  //  delete this->transform;
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
  
  this->density_estimator.reset();
  
  this->subsample_density_estimator.reset();
  
  this->the_worker.reset();
  
  this->proposal.reset();
  
  this->subsample_proposal.reset();
  
  //if (this->transform!=NULL)
  //  delete this->transform;
  
  LikelihoodEstimator::operator=(another);
  this->make_copy(another);
}

LikelihoodEstimator* DensityLikelihoodEstimator::duplicate() const
{
  return( new DensityLikelihoodEstimator(*this));
}

void DensityLikelihoodEstimator::make_copy(const DensityLikelihoodEstimator &another)
{
  if (another.density_estimator!=NULL)
    this->density_estimator.reset(another.density_estimator->duplicate());
  else
    this->density_estimator.reset();
  
  if (another.subsample_density_estimator!=NULL)
    this->subsample_density_estimator.reset(another.subsample_density_estimator->duplicate());
  else
    this->subsample_density_estimator.reset();
  
  if (another.the_worker!=NULL)
    this->the_worker.reset(another.the_worker->duplicate());
  else
    this->the_worker.reset();
  
  if (another.proposal!=NULL)
    this->proposal.reset(another.proposal->independent_proposal_kernel_duplicate());
  else
    this->proposal.reset();
  
  if (another.subsample_proposal!=NULL)
    this->subsample_proposal.reset(another.subsample_proposal->independent_proposal_kernel_duplicate());
  else
    this->subsample_proposal.reset();
  
  this->store_output = another.store_output;
  
  /*
   if (another.transform!=NULL)
   this->transform = another.transform->duplicate();
   else
   this->transform = NULL;
   */
  
  //this->simulate_distribution = another.simulate_distribution;
  
  //this->subsample_simulate_distribution = another.subsample_simulate_distribution;
  
  this->number_of_points = another.number_of_points;
  
  this->variables = another.variables;
}

// double DensityLikelihoodEstimator::estimate_log_likelihood(const List &inputs,
//                                                          const List &auxiliary_variables) const
// {
//   return this->func(inputs,this->observed_data);
// }

LikelihoodEstimatorOutput* DensityLikelihoodEstimator::initialise()
{
  return new DensityLikelihoodEstimatorOutput(this,
                                              this->density_estimator.get(),
                                              this->subsample_density_estimator.get(),
                                              this->store_output);
}

LikelihoodEstimatorOutput* DensityLikelihoodEstimator::initialise(const Parameters &parameters)
{
  return new DensityLikelihoodEstimatorOutput(this,
                                              this->density_estimator.get(),
                                              this->subsample_density_estimator.get(),
                                              this->store_output);
}

void DensityLikelihoodEstimator::setup()
{
  this->variables = this->proposal->independent_simulate(*this->rng).get_vector_variables();
  /*
   if (this->transform==NULL)
   this->variables = this->proposal->independent_simulate(*this->rng).get_vector_variables();
   else
   this->variables = this->transform->transform(this->proposal->independent_simulate(*this->rng)).get_vector_variables();
   */
}

void DensityLikelihoodEstimator::setup(const Parameters &parameters)
{
  /*
   if (this->transform==NULL)
   this->variables = this->proposal->independent_simulate(*this->rng,
   parameters).get_vector_variables();
   else
   this->variables = this->transform->transform(this->proposal->independent_simulate(*this->rng,
   parameters)).get_vector_variables();
   */
  this->variables = this->proposal->independent_simulate(*this->rng,
                                                         parameters).get_vector_variables();
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

void DensityLikelihoodEstimator::specific_change_data(Data* new_data)
{
}
}
