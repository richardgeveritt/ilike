#ifndef MOVEOUTPUT_H
#define MOVEOUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <iostream>
#include "particle.h"
#include "ilike_header.h"

namespace ilike
{
  /**
   * @file move_output.h
   * @brief Defines the Ensemble class.
   *
   * Provides ensemble functionality.
   *
   * @namespace ilike
   * @class Ensemble
   * @brief The ensemble class.
   */


class Ensemble;
class Particles;

class MoveOutput
{
  
public:
  
  /**
   * @brief Performs the moveoutput operation.
   */
  MoveOutput();
  /**
   * @brief Performs the ~moveoutput operation.
   */
  virtual ~MoveOutput();
  
  /**
   * @brief Performs the moveoutput operation.
   *
   * @param another The Ensemble instance to copy from.
   */
  MoveOutput(const MoveOutput &another);
  
  /**
   * @brief Assignment operator for Ensemble.
   *
   * @param another The Ensemble instance to copy from.
   */
  void operator=(const MoveOutput &another);
  /**
   * @brief Creates a deep copy of this Ensemble object.
   *
   * @return The result.
   */
  virtual MoveOutput* duplicate() const=0;
  
  /**
   * @brief Performs the back operation.
   *
   * @return The result.
   */
  virtual Particle& back()=0;
  /**
   * @brief Performs the back operation.
   *
   * @return The result.
   */
  virtual Particle back() const=0;
  
  /**
   * @brief Returns the vector of parameters.
   *
   * @return The result.
   */
  virtual std::vector<Parameters> get_vector_of_parameters() const=0;
  
  virtual void write_vector_points(const std::vector<std::string> &variables,
                                   std::ofstream &file_stream,
                                   std::shared_ptr<Transform> transform) const=0;
  virtual void write_any_points(const std::vector<std::string> &variables,
                                std::ofstream &file_stream) const=0;

  // Returns all vector points as a matrix (one row per MCMC step / particle).
  virtual arma::mat get_matrix_of_vector_points(const std::vector<std::string> &variables,
                                                std::shared_ptr<Transform> transform) const=0;
  
  virtual void write_factors(const std::string &directory_name,
                             const std::string &index) const=0;
  virtual void write_ensemble_factors(const std::string &directory_name,
                                      const std::string &index) const=0;
  
  /**
   * @brief Performs the length operation.
   *
   * @return The result.
   */
  virtual size_t length() const=0;
  
  /**
   * @brief Closes any open file streams.
   */
  virtual void close_ofstreams()=0;
  
  /**
   * @brief Returns the current algorithm parameters.
   *
   * @return The result.
   */
  virtual Parameters get_current_algorithm_parameters() const=0;
  
protected:
  
  /**
   * @brief Copies the state of another Ensemble into this object.
   *
   * @param another The Ensemble instance to copy from.
   */
  void make_copy(const MoveOutput &another);
  
  //std::vector<std::string> vector_variables;
  //std::vector<std::string> any_variables;
  
};
}

#endif
