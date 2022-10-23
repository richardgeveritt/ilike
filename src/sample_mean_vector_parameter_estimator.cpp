#include "sample_mean_vector_parameter_estimator.h"
#include "utils.h"
#include "move_output.h"

SampleMeanVectorParameterEstimator::SampleMeanVectorParameterEstimator()
  :VectorParameterEstimator()
{
}

SampleMeanVectorParameterEstimator::~SampleMeanVectorParameterEstimator()
{
  
}

//Copy constructor for the SampleMeanVectorParameterEstimator class.
SampleMeanVectorParameterEstimator::SampleMeanVectorParameterEstimator(const SampleMeanVectorParameterEstimator &another)
  :VectorParameterEstimator(another)
{
  this->make_copy(another);
}

void SampleMeanVectorParameterEstimator::operator=(const SampleMeanVectorParameterEstimator &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  VectorParameterEstimator::operator=(another);
  this->make_copy(another);
}

VectorParameterEstimator* SampleMeanVectorParameterEstimator::duplicate(void) const
{
  return( new SampleMeanVectorParameterEstimator(*this));
}

void SampleMeanVectorParameterEstimator::make_copy(const SampleMeanVectorParameterEstimator &another)
{
}

void SampleMeanVectorParameterEstimator::fit(const std::vector<Parameters> &points,
                                             arma::colvec normalised_log_weights)
{
  arma::mat matrix_points = vector_of_parameters_to_mat(points);
  arma::colvec wt = exp(normalised_log_weights);
  this->estimated = mean_wt(matrix_points,wt);
}

void SampleMeanVectorParameterEstimator::fit(const std::string &variable,
                                             const std::vector<Parameters> &points,
                                             arma::colvec normalised_log_weights)
{
  arma::mat matrix_points = vector_of_parameters_to_mat(variable,
                                                        points);
  arma::colvec wt = exp(normalised_log_weights);
  this->estimated = mean_wt(matrix_points,wt);
}

void SampleMeanVectorParameterEstimator::fit(const std::string &variable,
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
  this->estimated = mean_wt(matrix_points,wt);
}
