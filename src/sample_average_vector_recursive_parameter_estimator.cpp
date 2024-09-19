#include "sample_average_vector_recursive_parameter_estimator.h"
#include "utils.h"
#include "proposal_kernel.h"

namespace ilike
{
SampleAverageVectorRecursiveParameterEstimator::SampleAverageVectorRecursiveParameterEstimator()
:VectorRecursiveParameterEstimator()
{
  this->gain = equal_weight_gain;
}

SampleAverageVectorRecursiveParameterEstimator::~SampleAverageVectorRecursiveParameterEstimator()
{
  
}

//Copy constructor for the SampleAverageVectorRecursiveParameterEstimator class.
SampleAverageVectorRecursiveParameterEstimator::SampleAverageVectorRecursiveParameterEstimator(const SampleAverageVectorRecursiveParameterEstimator &another)
:VectorRecursiveParameterEstimator(another)
{
  this->make_copy(another);
}

void SampleAverageVectorRecursiveParameterEstimator::operator=(const SampleAverageVectorRecursiveParameterEstimator &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  VectorRecursiveParameterEstimator::operator=(another);
  this->make_copy(another);
}

RecursiveParameterEstimator* SampleAverageVectorRecursiveParameterEstimator::duplicate(void) const
{
  return( new SampleAverageVectorRecursiveParameterEstimator(*this));
}

VectorRecursiveParameterEstimator* SampleAverageVectorRecursiveParameterEstimator::vector_duplicate(void) const
{
  return( new SampleAverageVectorRecursiveParameterEstimator(*this));
}

void SampleAverageVectorRecursiveParameterEstimator::make_copy(const SampleAverageVectorRecursiveParameterEstimator &another)
{
  this->gain = another.gain;
}

void SampleAverageVectorRecursiveParameterEstimator::update(const std::string &variable_name,
                                                            const Particle &latest_particle,
                                                            size_t iteration_counter,
                                                            ProposalKernel* proposal)
{
  Rcpp::stop("SampleAverageVectorRecursiveParameterEstimator::update - need to sort out transform.");
  
  /*
   if (iteration_counter==0)
   {
   this->estimated = latest_particle.move_parameters->get_colvec(variable_name);
   }
   else
   {
   this->estimated = this->estimated + this->gain(iteration_counter+1)*(latest_particle.move_parameters->get_colvec(variable_name) - this->estimated);
   }
   */
}
}
