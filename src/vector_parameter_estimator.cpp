#include "vector_parameter_estimator.h"
#include "shared_point_weights.h"
#include "utils.h"
#include "move_output.h"

namespace ilike
{
VectorParameterEstimator::VectorParameterEstimator()
:ParameterEstimator()
{
}

VectorParameterEstimator::~VectorParameterEstimator()
{
  
}

//Copy constructor for the VectorParameterEstimator class.
VectorParameterEstimator::VectorParameterEstimator(const VectorParameterEstimator &another)
:ParameterEstimator(another)
{
  this->make_copy(another);
}

void VectorParameterEstimator::operator=(const VectorParameterEstimator &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  ParameterEstimator::operator=(another);
  this->make_copy(another);
}

void VectorParameterEstimator::make_copy(const VectorParameterEstimator &another)
{
  this->estimated = another.estimated;
}

void VectorParameterEstimator::fit(const std::vector<std::string> &variables,
                                   const std::vector<Parameters> &points,
                                   const arma::colvec &normalised_log_weights)
{
  arma::mat matrix_points = vector_of_parameters_to_mat(variables,
                                                        points);
  this->fit(matrix_points,exp(normalised_log_weights));
}

void VectorParameterEstimator::fit(const std::string &variable,
                                   const std::vector<Parameters> &points,
                                   const arma::colvec &normalised_log_weights)
{
  arma::mat matrix_points = vector_of_parameters_to_mat(variable,
                                                        points);
  this->fit(matrix_points,exp(normalised_log_weights));
}

void VectorParameterEstimator::fit(const std::string &variable,
                                   const std::vector<MoveOutput*> &points,
                                   const arma::colvec &normalised_log_weights)
{
  // each MoveOutput is one particle, which may hold several points: its
  // weight is shared between them (shared_point_weights.h)
  std::vector< std::vector<Parameters> > points_per_particle;
  points_per_particle.reserve(points.size());
  for (auto i=points.begin(); i!=points.end(); ++i)
    points_per_particle.push_back((*i)->get_vector_of_parameters());
  std::vector<Parameters> all_points;
  arma::colvec all_normalised_log_weights;
  shared_point_weights(points_per_particle,
                       normalised_log_weights,
                       all_points,
                       all_normalised_log_weights);
  arma::mat matrix_points = vector_of_parameters_to_mat(variable,
                                                        all_points);
  this->fit(matrix_points,exp(all_normalised_log_weights));
}
}
