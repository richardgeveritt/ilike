#include "sequential_density_likelihood_estimator_worker.h"
#include "likelihood_estimator_output.h"
#include "density_likelihood_estimator.h"
#include "utils.h"

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
SequentialDensityLikelihoodEstimatorWorker::~SequentialDensityLikelihoodEstimatorWorker(void)
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

void SequentialDensityLikelihoodEstimatorWorker::specific_simulate()
{
  RandomNumberGenerator local_rng(*this->get_rng());
  local_rng.seed(this->get_seed(),this->get_number_of_points());
  
  this->points.clear();
  this->points.reserve(this->get_number_of_points());
  
  for (size_t i = 0; i < this->get_number_of_points(); ++i)
  {
    this->points.push_back(this->the_dle->simulate_distribution(local_rng,Parameters()));
  }
}

void SequentialDensityLikelihoodEstimatorWorker::specific_simulate(const Parameters &conditioned_on_parameters)
{
  RandomNumberGenerator local_rng(*this->get_rng());
  local_rng.seed(this->get_seed(),this->get_number_of_points());
  
  this->points.clear();
  this->points.reserve(this->get_number_of_points());
  
  for (size_t i = 0; i < this->get_number_of_points(); ++i)
  {
    this->points.push_back(this->the_dle->simulate_distribution(local_rng,
                                                                conditioned_on_parameters));
  }
}

void SequentialDensityLikelihoodEstimatorWorker::subsample_specific_simulate(const Parameters &conditioned_on_parameters)
{
  RandomNumberGenerator local_rng(*this->get_rng());
  local_rng.seed(this->get_seed(),this->get_number_of_points());
  
  this->points.clear();
  this->points.reserve(this->get_number_of_points());
  
  for (size_t i = 0; i < this->get_number_of_points(); ++i)
  {
    this->points.push_back(this->the_dle->subsample_simulate_distribution(local_rng,
                                                                conditioned_on_parameters));
  }
}
