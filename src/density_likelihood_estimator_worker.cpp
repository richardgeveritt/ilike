#include "density_likelihood_estimator_worker.h"
#include "likelihood_estimator.h"
#include "likelihood_estimator_output.h"
#include "density_likelihood_estimator.h"

namespace ilike
{
DensityLikelihoodEstimatorWorker::DensityLikelihoodEstimatorWorker()
{
}

// Constructor for IS where prior is proposal.
DensityLikelihoodEstimatorWorker::DensityLikelihoodEstimatorWorker(DensityLikelihoodEstimator* the_dle_in)
{
  this->the_dle = the_dle_in;
}

DensityLikelihoodEstimatorWorker::DensityLikelihoodEstimatorWorker(const DensityLikelihoodEstimatorWorker &another)
{
  this->make_copy(another);
}

DensityLikelihoodEstimatorWorker::~DensityLikelihoodEstimatorWorker()
{
}

void DensityLikelihoodEstimatorWorker::operator=(const DensityLikelihoodEstimatorWorker &another)
{
  if(this == &another)
    return;
  
  this->make_copy(another);
}

void DensityLikelihoodEstimatorWorker::make_copy(const DensityLikelihoodEstimatorWorker &another)
{
  this->the_dle = another.the_dle;

}

void DensityLikelihoodEstimatorWorker::simulate(DensityLikelihoodEstimatorOutput* output)
{
  this->specific_simulate(output);
  this->set_seed(this->get_seed() + this->get_number_of_points());
}

void DensityLikelihoodEstimatorWorker::simulate(DensityLikelihoodEstimatorOutput* output,
                                                const Parameters &conditioned_on_parameters)
{
  this->specific_simulate(output,
                          conditioned_on_parameters);
  this->set_seed(this->get_seed() + this->get_number_of_points());
}

void DensityLikelihoodEstimatorWorker::subsample_simulate(DensityLikelihoodEstimatorOutput* output)
{
  this->subsample_specific_simulate(output);
  this->set_seed(this->get_seed() + this->get_number_of_points());
}

void DensityLikelihoodEstimatorWorker::subsample_simulate(DensityLikelihoodEstimatorOutput* output,
                                                          const Parameters &conditioned_on_parameters)
{
  this->subsample_specific_simulate(output,
                                    conditioned_on_parameters);
  this->set_seed(this->get_seed() + this->get_number_of_points());
}

size_t DensityLikelihoodEstimatorWorker::get_number_of_points() const
{
  return this->the_dle->number_of_points;
}

RandomNumberGenerator* DensityLikelihoodEstimatorWorker::get_rng()
{
  return this->the_dle->rng;
}

size_t DensityLikelihoodEstimatorWorker::get_seed() const
{
  return *this->the_dle->seed;
}

void DensityLikelihoodEstimatorWorker::set_seed(size_t seed_in)
{
  *this->the_dle->seed = seed_in;
}
}
