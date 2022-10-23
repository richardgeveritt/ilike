#include "matrix_recursive_parameter_estimator.h"

MatrixRecursiveParameterEstimator::MatrixRecursiveParameterEstimator()
  :RecursiveParameterEstimator()
{
}

MatrixRecursiveParameterEstimator::~MatrixRecursiveParameterEstimator()
{
  
}

//Copy constructor for the MatrixRecursiveParameterEstimator class.
MatrixRecursiveParameterEstimator::MatrixRecursiveParameterEstimator(const MatrixRecursiveParameterEstimator &another)
  :RecursiveParameterEstimator(another)
{
  this->make_copy(another);
}

void MatrixRecursiveParameterEstimator::operator=(const MatrixRecursiveParameterEstimator &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  RecursiveParameterEstimator::operator=(another);
  this->make_copy(another);
}

void MatrixRecursiveParameterEstimator::make_copy(const MatrixRecursiveParameterEstimator &another)
{
  this->estimated = another.estimated;
}
