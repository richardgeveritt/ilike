#ifndef PDMPMCMCOUTPUT_H
#define PDMPMCMCOUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "move_output.h"

namespace ilike
{
  /**
   * @file pdmp_mcmc_output.h
   * @brief Defines the PDMPMCMCOutput class.
   *
   * Stores and manages the output produced by PDMPMCMC. Holds results such as log-likelihood estimates, samples, or diagnostics returned after running the associated algorithm.
   *
   * @namespace ilike
   * @class PDMPMCMCOutput
   * @brief A pdmpmcmc output derived from MoveOutput.
   */


class PDMPMCMCOutput : public MoveOutput
{
  
public:
  
  /**
   * @brief Default constructor for PDMPMCMCOutput.
   */
  PDMPMCMCOutput();
  
  /**
   * @brief Destructor for PDMPMCMCOutput.
   */
  virtual ~PDMPMCMCOutput();
  
  /**
   * @brief Copy constructor for PDMPMCMCOutput.
   *
   * @param another The PDMPMCMCOutput instance to copy from.
   */
  PDMPMCMCOutput(const PDMPMCMCOutput &another);
  
  /**
   * @brief Assignment operator for PDMPMCMCOutput.
   *
   * @param another The PDMPMCMCOutput instance to copy from.
   */
  void operator=(const PDMPMCMCOutput &another);
  /**
   * @brief Creates a deep copy of this PDMPMCMCOutput object.
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
                           std::shared_ptr<Transform> transform) const;
  void write_any_points(const std::vector<std::string> &variables,
                        std::ofstream &file_stream) const;

  arma::mat get_matrix_of_vector_points(const std::vector<std::string> &variables,
                                        std::shared_ptr<Transform> transform) const;
  
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
  
  /** @brief The dummy. */
  Particle dummy;
  
  /**
   * @brief Copies the state of another PDMPMCMCOutput into this object.
   *
   * @param another The PDMPMCMCOutput instance to copy from.
   */
  void make_copy(const PDMPMCMCOutput &another);
  
};
}

#endif
