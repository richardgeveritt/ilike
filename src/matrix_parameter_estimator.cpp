#include "matrix_parameter_estimator.h"

MatrixParameterEstimator::MatrixParameterEstimator()
  :ParameterEstimator()
{
}

MatrixParameterEstimator::~MatrixParameterEstimator()
{
  
}

//Copy constructor for the MatrixParameterEstimator class.
MatrixParameterEstimator::MatrixParameterEstimator(const MatrixParameterEstimator &another)
  :ParameterEstimator(another)
{
  this->make_copy(another);
}

void MatrixParameterEstimator::operator=(const MatrixParameterEstimator &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  ParameterEstimator::operator=(another);
  this->make_copy(another);
}

void MatrixParameterEstimator::make_copy(const MatrixParameterEstimator &another)
{
  this->estimated = another.estimated;
}
