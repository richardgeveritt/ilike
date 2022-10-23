#include "sample_covariance_matrix_parameter_estimator.h"
#include "utils.h"
#include "move_output.h"

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

void SampleCovarianceMatrixParameterEstimator::fit(const std::vector<Parameters> &points,
                                                   arma::colvec normalised_log_weights)
{
  arma::mat matrix_points = vector_of_parameters_to_mat(points);
  arma::colvec wt = exp(normalised_log_weights);
  this->estimated = cov_wt(matrix_points,wt);
}

void SampleCovarianceMatrixParameterEstimator::fit(const std::string &variable,
                                                   const std::vector<Parameters> &points,
                                                   arma::colvec normalised_log_weights)
{
  arma::mat matrix_points = vector_of_parameters_to_mat(variable,
                                                        points);
  arma::colvec wt = exp(normalised_log_weights);
  this->estimated = cov_wt(matrix_points,wt);
}

void SampleCovarianceMatrixParameterEstimator::fit(const std::string &variable,
                                             const std::vector<MoveOutput*> &points,
                                             arma::colvec normalised_log_weights)
{
  std::vector<Parameters> all_points;
  std::vector<Parameters> current_parameters = (*points.begin())->get_vector_of_parameters();
  // assume that all MoveOutput* give the same number of points.
  all_points.reserve(points.size()*current_parameters.size());
  arma::colvec all_normalised_log_weights(points.size()*current_parameters.size());
  
  size_t total_counter = 0;
  size_t outer_counter = 0;
  for (auto i=points.begin(); i!=points.end(); ++i)
  {
    std::vector<Parameters> current_parameters = (*i)->get_vector_of_parameters();
    for (auto j=current_parameters.begin(); j!=current_parameters.end(); ++j)
    {
      all_points.push_back(*j);
      all_normalised_log_weights[total_counter] = normalised_log_weights[outer_counter];
    }
    outer_counter = outer_counter + 1;
  }
  arma::mat matrix_points = vector_of_parameters_to_mat(variable,
                                                        all_points);
  arma::colvec wt = exp(all_normalised_log_weights);
  this->estimated = cov_wt(matrix_points,wt);
}
