#include "matrix_parameter_estimator.h"
#include "utils.h"
#include "move_output.h"

namespace ilike
{
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

void MatrixParameterEstimator::fit(const std::vector<std::string> &variables,
                                   const std::vector<Parameters> &points,
                                   const arma::colvec &normalised_log_weights)
{
  arma::mat matrix_points = vector_of_parameters_to_mat(variables,
                                                        points);
  this->fit(matrix_points,exp(normalised_log_weights));
}

void MatrixParameterEstimator::fit(const std::string &variable,
                                   const std::vector<Parameters> &points,
                                   const arma::colvec &normalised_log_weights)
{
  arma::mat matrix_points = vector_of_parameters_to_mat(variable,
                                                        points);
  this->fit(matrix_points,exp(normalised_log_weights));
}

void MatrixParameterEstimator::fit(const std::string &variable,
                                   const std::vector<MoveOutput*> &points,
                                   const arma::colvec &normalised_log_weights)
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
  this->fit(matrix_points,exp(normalised_log_weights));
}
}
