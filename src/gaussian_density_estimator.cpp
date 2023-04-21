#include "gaussian_density_estimator.h"
#include "sample_mean_vector_parameter_estimator.h"
#include "sample_covariance_matrix_parameter_estimator.h"
#include "distributions.h"
#include "gaussian_density_estimator_output.h"

GaussianDensityEstimator::GaussianDensityEstimator()
  :DensityEstimator()
{
  this->unbiased = true;
}

GaussianDensityEstimator::GaussianDensityEstimator(const std::vector<std::string> &variables_in)
:DensityEstimator(variables_in)
{
  //this->mean_estimator = new SampleMeanVectorParameterEstimator();
  //this->covariance_estimator = new SampleCovarianceMatrixParameterEstimator();
}

GaussianDensityEstimator::GaussianDensityEstimator(const std::vector<std::string> &variables_in,
                                                   bool unbiased_in)
  :DensityEstimator(variables_in)
{
  this->unbiased = unbiased_in;
  //this->variables = variables_in;
  
  //this->mean_estimator = new SampleMeanVectorParameterEstimator();
  //this->covariance_estimator = new SampleCovarianceMatrixParameterEstimator();
}

GaussianDensityEstimator::~GaussianDensityEstimator()
{
  //if (this->mean_estimator==NULL)
  //  delete this->mean_estimator;
  
  //if (this->covariance_estimator==NULL)
  //  delete this->covariance_estimator;
}

//Copy constructor for the GaussianDensityEstimator class.
GaussianDensityEstimator::GaussianDensityEstimator(const GaussianDensityEstimator &another)
  :DensityEstimator(another)
{
  this->make_copy(another);
}

void GaussianDensityEstimator::operator=(const GaussianDensityEstimator &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  //if (this->mean_estimator==NULL)
  //  delete this->mean_estimator;
  
  //if (this->covariance_estimator==NULL)
  //  delete this->covariance_estimator;
  
  DensityEstimator::operator=(another);
  this->make_copy(another);
}

DensityEstimator* GaussianDensityEstimator::duplicate() const
{
  return( new GaussianDensityEstimator(*this));
}

void GaussianDensityEstimator::make_copy(const GaussianDensityEstimator &another)
{
  /*
  if (another.mean_estimator!=NULL)
    this->mean_estimator = another.mean_estimator->duplicate();
  else
    this->mean_estimator = NULL;
  
  if (another.covariance_estimator!=NULL)
    this->covariance_estimator = another.covariance_estimator->duplicate();
  else
    this->covariance_estimator = NULL;
  */
  this->unbiased = another.unbiased;
  
}

DensityEstimatorOutput* GaussianDensityEstimator::initialise()
{
  return new GaussianDensityEstimatorOutput(this);
}

bool GaussianDensityEstimator::get_unbiased() const
{
  return this->unbiased;
}

/*
 void GaussianDensityEstimator::fit(const std::vector<Parameters> &points,
 arma::colvec normalised_log_weights)
 {
 this->n = points.size();
 //this->mean_estimator->fit(this->variables, points, normalised_log_weights);
 //this->covariance_estimator->fit(this->variables, points, normalised_log_weights);
 }
 
 double GaussianDensityEstimator::evaluate(const Data &point) const
 {
 size_t d = this->mean_estimator->estimated.n_rows;
 if ( (unbiased==true) && (this->n>d+3) )
 {
 return(dmvnorm_estimated_params(point.get_colvec(this->variables),
 this->mean_estimator->estimated,
 this->covariance_estimator->estimated,
 this->n));
 }
 else
 {
 return(dmvnorm(point.get_colvec(this->variables),
 this->mean_estimator->estimated,
 this->covariance_estimator->estimated));
 }
 }
*/
