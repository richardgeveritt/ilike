#include "double_scale_recursive_parameter_estimator.h"
#include "double_recursive_parameter_estimator.h"
#include "utils.h"

DoubleScaleRecursiveParameterEstimator::DoubleScaleRecursiveParameterEstimator()
  :ScaleRecursiveParameterEstimator()
{
}

DoubleScaleRecursiveParameterEstimator::~DoubleScaleRecursiveParameterEstimator()
{
  if (this->recursive_estimator!=NULL)
    delete this->recursive_estimator;
}

//Copy constructor for the DoubleScaleRecursiveParameterEstimator class.
DoubleScaleRecursiveParameterEstimator::DoubleScaleRecursiveParameterEstimator(const DoubleScaleRecursiveParameterEstimator &another)
  :ScaleRecursiveParameterEstimator(another)
{
  this->make_copy(another);
}

void DoubleScaleRecursiveParameterEstimator::operator=(const DoubleScaleRecursiveParameterEstimator &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  if (this->recursive_estimator!=NULL)
    delete this->recursive_estimator;
  
  ScaleRecursiveParameterEstimator::operator=(another);
  this->make_copy(another);
}

RecursiveParameterEstimator* DoubleScaleRecursiveParameterEstimator::duplicate() const
{
  return( new DoubleScaleRecursiveParameterEstimator(*this));
}

ScaleRecursiveParameterEstimator* DoubleScaleRecursiveParameterEstimator::scale_duplicate() const
{
  return( new DoubleScaleRecursiveParameterEstimator(*this));
}

void DoubleScaleRecursiveParameterEstimator::make_copy(const DoubleScaleRecursiveParameterEstimator &another)
{
  if (another.recursive_estimator!=NULL)
    this->recursive_estimator = another.recursive_estimator->double_duplicate();
  else
    this->recursive_estimator = NULL;
}

void DoubleScaleRecursiveParameterEstimator::update(const std::string &variable_name,
                                                    const Particle &latest_particle,
                                                    size_t iteration_counter,
                                                    ProposalKernel* proposal)
{
  this->recursive_estimator->update(variable_name,
                                    latest_particle,
                                    iteration_counter,
                                    proposal);
  this->estimated.get_constant() = exp(this->recursive_estimator->estimated);
}
