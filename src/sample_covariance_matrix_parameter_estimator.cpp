#include "sample_covariance_matrix_parameter_estimator.h"
#include "utils.h"

namespace ilike
{
SampleCovarianceMatrixParameterEstimator::SampleCovarianceMatrixParameterEstimator()
:MatrixParameterEstimator()
{
}

SampleCovarianceMatrixParameterEstimator::~SampleCovarianceMatrixParameterEstimator()
{
  
}

//Copy constructor for the SampleCovarianceMatrixParameterEstimator class.
SampleCovarianceMatrixParameterEstimator::SampleCovarianceMatrixParameterEstimator(const SampleCovarianceMatrixParameterEstimator &another)
:MatrixParameterEstimator(another)
{
  this->make_copy(another);
}

void SampleCovarianceMatrixParameterEstimator::operator=(const SampleCovarianceMatrixParameterEstimator &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  MatrixParameterEstimator::operator=(another);
  this->make_copy(another);
}

MatrixParameterEstimator* SampleCovarianceMatrixParameterEstimator::duplicate(void) const
{
  return( new SampleCovarianceMatrixParameterEstimator(*this));
}

void SampleCovarianceMatrixParameterEstimator::make_copy(const SampleCovarianceMatrixParameterEstimator &another)
{
}

void SampleCovarianceMatrixParameterEstimator::fit(const arma::mat &matrix_points,
                                                   const arma::colvec &wt)
{
  this->estimated = cov_wt(matrix_points,wt);
  if (!this->estimated.is_sympd())
    this->estimated = this->estimated + 0.0000000000001*arma::eye(arma::size(this->estimated));
}
}
