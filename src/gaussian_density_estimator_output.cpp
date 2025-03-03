#include "gaussian_density_estimator_output.h"
#include "gaussian_density_estimator.h"
#include "sample_mean_vector_parameter_estimator.h"
#include "sample_covariance_matrix_parameter_estimator.h"
#include "distributions.h"
#include "utils.h"

namespace ilike
{
GaussianDensityEstimatorOutput::GaussianDensityEstimatorOutput()
:DensityEstimatorOutput()
{
  this->mean_estimator = NULL;
  this->covariance_estimator = NULL;
}

GaussianDensityEstimatorOutput::GaussianDensityEstimatorOutput(GaussianDensityEstimator* estimator_in)
:DensityEstimatorOutput()
{
  this->estimator = estimator_in;
  this->mean_estimator = new SampleMeanVectorParameterEstimator();
  this->covariance_estimator = new SampleCovarianceMatrixParameterEstimator();
}
GaussianDensityEstimatorOutput::~GaussianDensityEstimatorOutput()
{
  if (this->mean_estimator==NULL)
    delete this->mean_estimator;
  
  if (this->covariance_estimator==NULL)
    delete this->covariance_estimator;
}

//Copy constructor for the GaussianDensityEstimatorOutput class.
GaussianDensityEstimatorOutput::GaussianDensityEstimatorOutput(const GaussianDensityEstimatorOutput &another)
:DensityEstimatorOutput(another)
{
  this->make_copy(another);
}

void GaussianDensityEstimatorOutput::operator=(const GaussianDensityEstimatorOutput &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  if (this->mean_estimator==NULL)
    delete this->mean_estimator;
  
  if (this->covariance_estimator==NULL)
    delete this->covariance_estimator;
  
  DensityEstimatorOutput::operator=(another);
  this->make_copy(another);
}

DensityEstimatorOutput* GaussianDensityEstimatorOutput::duplicate() const
{
  return( new GaussianDensityEstimatorOutput(*this));
}

void GaussianDensityEstimatorOutput::make_copy(const GaussianDensityEstimatorOutput &another)
{
  if (another.mean_estimator!=NULL)
    this->mean_estimator = another.mean_estimator->duplicate();
  else
    this->mean_estimator = NULL;
  
  if (another.covariance_estimator!=NULL)
    this->covariance_estimator = another.covariance_estimator->duplicate();
  else
    this->covariance_estimator = NULL;
  
  //this->unbiased = another.unbiased;
  
}

void GaussianDensityEstimatorOutput::fit(const std::vector<Parameters> &points,
                                         const arma::colvec &normalised_log_weights)
{
  this->n = points.size();
  arma::mat matrix_points = vector_of_parameters_to_mat(this->estimator->get_variables(),
                                                        points);
  this->mean_estimator->fit(matrix_points, exp(normalised_log_weights));
  this->covariance_estimator->fit(matrix_points, exp(normalised_log_weights));
}

double GaussianDensityEstimatorOutput::evaluate(const Data &point) const
{
  size_t d = this->mean_estimator->estimated.n_rows;
  if ( (this->estimator->get_unbiased()==true) && (this->n>d+3) )
  {
    return(dmvnorm_estimated_params(point.get_colvec(this->estimator->get_variables()),
                                    this->mean_estimator->estimated,
                                    this->covariance_estimator->estimated,
                                    this->n));
  }
  else
  {

    return(dmvnorm(point.get_colvec(this->estimator->get_variables()),
                   this->mean_estimator->estimated,
                   this->covariance_estimator->estimated));
  }
}
}
