#ifndef SINGLEPOINTMOVEOUTPUT_H
#define SINGLEPOINTMOVEOUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "move_output.h"
#include "parameters.h"
#include "particle.h"
#include "ilike_header.h"

namespace ilike
{
  /**
   * @file single_point_move_output.h
   * @brief Defines the Factors class.
   *
   * Provides factors functionality.
   *
   * @namespace ilike
   * @class Factors
   * @brief The factors class.
   */


class Factors;

class SinglePointMoveOutput : public MoveOutput
{
  
public:
  
  SinglePointMoveOutput(const Parameters &parameter_in,
                        const Factors* factors_in,
                        const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
                        const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in);
  SinglePointMoveOutput(const Parameters &parameter_in,
                        const EnsembleFactors* factors_in);
  
  SinglePointMoveOutput(Parameters &&parameter_in,
                        const Factors* factors_in,
                        const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
                        const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in);
  SinglePointMoveOutput(Parameters &&parameter_in,
                        const EnsembleFactors* factors_in);
  
  /**
   * @brief Performs the singlepointmoveoutput operation.
   *
   * @param particle_in The particle.
   */
  SinglePointMoveOutput(const Particle &particle_in);
  /**
   * @brief Performs the singlepointmoveoutput operation.
   *
   * @param particle_in The particle.
   */
  SinglePointMoveOutput(Particle &&particle_in);
  /**
   * @brief Performs the singlepointmoveoutput operation.
   */
  SinglePointMoveOutput();
  
  /**
   * @brief Performs the ~singlepointmoveoutput operation.
   */
  virtual ~SinglePointMoveOutput();
  
  /**
   * @brief Performs the singlepointmoveoutput operation.
   *
   * @param another The Factors instance to copy from.
   */
  SinglePointMoveOutput(const SinglePointMoveOutput &another);
  
  /**
   * @brief Assignment operator for Factors.
   *
   * @param another The Factors instance to copy from.
   */
  void operator=(const SinglePointMoveOutput &another);
  /**
   * @brief Creates a deep copy of this Factors object.
   *
   * @return The result.
   */
  MoveOutput* duplicate() const;
  
  /**
   * @brief Performs the back operation.
   *
   * @return The result.
   */
  Particle& back();
  /**
   * @brief Performs the back operation.
   *
   * @return The result.
   */
  Particle back() const;
  
  /**
   * @brief Returns the vector of parameters.
   *
   * @return The result.
   */
  std::vector<Parameters> get_vector_of_parameters() const;
  
  void write_vector_points(const std::vector<std::string> &variables,
                           std::ofstream &file_stream,
                           std::shared_ptr<Transform> inverse_transform) const;
  void write_any_points(const std::vector<std::string> &variables,
                        std::ofstream &file_stream) const;
  
  void write_factors(const std::string &directory_name,
                     const std::string &index) const;
  void write_ensemble_factors(const std::string &directory_name,
                              const std::string &index) const;
  
  /**
   * @brief Performs the length operation.
   *
   * @return The result.
   */
  size_t length() const;
  
  /**
   * @brief Closes any open file streams.
   */
  void close_ofstreams();
  
  /**
   * @brief Returns the current algorithm parameters.
   *
   * @return The result.
   */
  Parameters get_current_algorithm_parameters() const;
  
protected:
  
  /**
   * @brief Copies the state of another Factors into this object.
   *
   * @param another The Factors instance to copy from.
   */
  void make_copy(const SinglePointMoveOutput &another);
  
  /** @brief The output. */
  Particle output;
  
  /** @brief The algorithm parameters. */
  Parameters algorithm_parameters;
  
};
}

#endif
