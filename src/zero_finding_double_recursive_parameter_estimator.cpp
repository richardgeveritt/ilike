#include "zero_finding_double_recursive_parameter_estimator.h"
#include "utils.h"

ZeroFindingDoubleRecursiveParameterEstimator::ZeroFindingDoubleRecursiveParameterEstimator()
  :DoubleRecursiveParameterEstimator()
{
}

ZeroFindingDoubleRecursiveParameterEstimator::ZeroFindingDoubleRecursiveParameterEstimator(double initial_value,
                                                                                           double target_score_in)
:DoubleRecursiveParameterEstimator(initial_value)
{
  this->target_score = target_score_in;
}

ZeroFindingDoubleRecursiveParameterEstimator::~ZeroFindingDoubleRecursiveParameterEstimator()
{
  
}

//Copy constructor for the ZeroFindingDoubleRecursiveParameterEstimator class.
ZeroFindingDoubleRecursiveParameterEstimator::ZeroFindingDoubleRecursiveParameterEstimator(const ZeroFindingDoubleRecursiveParameterEstimator &another)
  :DoubleRecursiveParameterEstimator(another)
{
  this->make_copy(another);
}

void ZeroFindingDoubleRecursiveParameterEstimator::operator=(const ZeroFindingDoubleRecursiveParameterEstimator &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  DoubleRecursiveParameterEstimator::operator=(another);
  this->make_copy(another);
}

RecursiveParameterEstimator* ZeroFindingDoubleRecursiveParameterEstimator::duplicate(void) const
{
  return( new ZeroFindingDoubleRecursiveParameterEstimator(*this));
}

DoubleRecursiveParameterEstimator* ZeroFindingDoubleRecursiveParameterEstimator::double_duplicate(void) const
{
  return( new ZeroFindingDoubleRecursiveParameterEstimator(*this));
}

void ZeroFindingDoubleRecursiveParameterEstimator::make_copy(const ZeroFindingDoubleRecursiveParameterEstimator &another)
{
  this->gain = another.gain;
}

void ZeroFindingDoubleRecursiveParameterEstimator::update(const std::string &variable_name,
                                                          const Particle &latest_particle,
                                                          size_t iteration_counter,
                                                          ProposalKernel* proposal)
{
  // only uses acceptance at the moment
  this->estimated = this->estimated + this->gain(iteration_counter+1)*(double(latest_particle.accepted_outputs.find(proposal)->second) - this->target_score);
}
