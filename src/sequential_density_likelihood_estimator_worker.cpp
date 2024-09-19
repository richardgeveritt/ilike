#include "sequential_density_likelihood_estimator_worker.h"
#include "likelihood_estimator_output.h"
#include "density_likelihood_estimator.h"
#include "utils.h"
#include "independent_proposal_kernel.h"
#include "transform.h"
#include "density_likelihood_estimator_output.h"

namespace ilike
{
//Default constructor.
SequentialDensityLikelihoodEstimatorWorker::SequentialDensityLikelihoodEstimatorWorker()
:DensityLikelihoodEstimatorWorker()
{
}

SequentialDensityLikelihoodEstimatorWorker::SequentialDensityLikelihoodEstimatorWorker(DensityLikelihoodEstimator* the_dle_in)
:DensityLikelihoodEstimatorWorker(the_dle_in)
{
}

//Copy constructor for the SequentialDensityLikelihoodEstimatorWorker class.
SequentialDensityLikelihoodEstimatorWorker::SequentialDensityLikelihoodEstimatorWorker(const SequentialDensityLikelihoodEstimatorWorker &another)
:DensityLikelihoodEstimatorWorker(another)
{
  this->make_copy(another);
}

//Destructor for the SequentialDensityLikelihoodEstimatorWorker class.
SequentialDensityLikelihoodEstimatorWorker::~SequentialDensityLikelihoodEstimatorWorker()
{
}

void SequentialDensityLikelihoodEstimatorWorker::operator=(const SequentialDensityLikelihoodEstimatorWorker &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  DensityLikelihoodEstimatorWorker::operator=(another);
  this->make_copy(another);
}

DensityLikelihoodEstimatorWorker* SequentialDensityLikelihoodEstimatorWorker::duplicate() const
{
  return( new SequentialDensityLikelihoodEstimatorWorker(*this));
}

void SequentialDensityLikelihoodEstimatorWorker::make_copy(const SequentialDensityLikelihoodEstimatorWorker &another)
{
  this->points = another.points;
}

std::vector<Parameters> SequentialDensityLikelihoodEstimatorWorker::get_points() const
{
  return this->points;
}

void SequentialDensityLikelihoodEstimatorWorker::specific_simulate(DensityLikelihoodEstimatorOutput* output)
{
  RandomNumberGenerator local_rng(*this->get_rng());
  local_rng.seed(this->get_seed(),this->get_number_of_points());
  
  this->points.clear();
  this->points.reserve(this->get_number_of_points());
  
  for (size_t i = 0; i < this->get_number_of_points(); ++i)
  {
    this->points.push_back(this->the_dle->proposal->independent_simulate(local_rng));
  }
  
  output->fit(this->points);
  
  /*
   if (this->the_dle->transform==NULL)
   {
   for (size_t i = 0; i < this->get_number_of_points(); ++i)
   {
   this->points.push_back(this->the_dle->proposal->independent_simulate(local_rng));
   }
   }
   else
   {
   for (size_t i = 0; i < this->get_number_of_points(); ++i)
   {
   this->points.push_back(this->the_dle->transform->transform(this->the_dle->proposal->independent_simulate(local_rng)));
   }
   }
   */
}

void SequentialDensityLikelihoodEstimatorWorker::specific_simulate(DensityLikelihoodEstimatorOutput* output,
                                                                   const Parameters &conditioned_on_parameters)
{
  RandomNumberGenerator local_rng(*this->get_rng());
  local_rng.seed(this->get_seed(),this->get_number_of_points());
  
  this->points.clear();
  this->points.reserve(this->get_number_of_points());
  
  for (size_t i = 0; i < this->get_number_of_points(); ++i)
  {
    this->points.push_back(this->the_dle->proposal->independent_simulate(local_rng,
                                                                         conditioned_on_parameters));
  }
  
  output->fit(this->points);
  
  /*
   if (this->the_dle->transform==NULL)
   {
   for (size_t i = 0; i < this->get_number_of_points(); ++i)
   {
   this->points.push_back(this->the_dle->proposal->independent_simulate(local_rng,
   conditioned_on_parameters));
   }
   }
   else
   {
   for (size_t i = 0; i < this->get_number_of_points(); ++i)
   {
   this->points.push_back(this->the_dle->transform->transform(this->the_dle->proposal->independent_simulate(local_rng,
   conditioned_on_parameters)));
   }
   }
   */
}

void SequentialDensityLikelihoodEstimatorWorker::subsample_specific_simulate(DensityLikelihoodEstimatorOutput* output)
{
  RandomNumberGenerator local_rng(*this->get_rng());
  local_rng.seed(this->get_seed(),this->get_number_of_points());
  
  this->points.clear();
  this->points.reserve(this->get_number_of_points());
  
  for (size_t i = 0; i < this->get_number_of_points(); ++i)
  {
    this->points.push_back(this->the_dle->proposal->subsample_independent_simulate(local_rng));
  }
  
  output->subsample_fit(this->points);
  
  /*
   if (this->the_dle->transform==NULL)
   {
   
   }
   else
   {
   for (size_t i = 0; i < this->get_number_of_points(); ++i)
   {
   this->points.push_back(this->the_dle->transform->transform(this->the_dle->proposal->subsample_independent_simulate(local_rng)));
   }
   }
   */
}

void SequentialDensityLikelihoodEstimatorWorker::subsample_specific_simulate(DensityLikelihoodEstimatorOutput* output,
                                                                             const Parameters &conditioned_on_parameters)
{
  RandomNumberGenerator local_rng(*this->get_rng());
  local_rng.seed(this->get_seed(),this->get_number_of_points());
  
  this->points.clear();
  this->points.reserve(this->get_number_of_points());
  
  for (size_t i = 0; i < this->get_number_of_points(); ++i)
  {
    this->points.push_back(this->the_dle->proposal->subsample_independent_simulate(local_rng,
                                                                                   conditioned_on_parameters));
  }
  
  output->subsample_fit(this->points);
  
  /*
   if (this->the_dle->transform==NULL)
   {
   
   }
   else
   {
   for (size_t i = 0; i < this->get_number_of_points(); ++i)
   {
   this->points.push_back(this->the_dle->transform->transform(this->the_dle->proposal->subsample_independent_simulate(local_rng,
   conditioned_on_parameters)));
   }
   }
   */
}
}
