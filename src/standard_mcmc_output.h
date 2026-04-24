#ifndef STANDARDMCMCOUTPUT_H
#define STANDARDMCMCOUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "move_output.h"
#include "parameters.h"
#include "ilike_header.h"

namespace ilike
{
  /**
   * @file standard_mcmc_output.h
   * @brief Defines the MCMCTermination class.
   *
   * Implements an MCMC termination criterion based on default. The sampler queries this criterion after each step.
   *
   * @namespace ilike
   * @class MCMCTermination
   * @brief The mcmc termination class.
   */


class MCMCTermination;
class MCMC;

class StandardMCMCOutput : public MoveOutput
{
  
public:
  
  //StandardMCMCOutput(const Parameters &parameter_in);
  /**
   * @brief Performs the standardmcmcoutput operation.
   */
  StandardMCMCOutput();
  
  /**
   * @brief Performs the standardmcmcoutput operation.
   *
   * @param termination_in The termination.
   */
  StandardMCMCOutput(MCMCTermination* termination_in);
  
  /**
   * @brief Performs the ~standardmcmcoutput operation.
   */
  virtual ~StandardMCMCOutput();
  
  /**
   * @brief Performs the standardmcmcoutput operation.
   *
   * @param another The MCMCTermination instance to copy from.
   */
  StandardMCMCOutput(const StandardMCMCOutput &another);
  
  /**
   * @brief Assignment operator for MCMCTermination.
   *
   * @param another The MCMCTermination instance to copy from.
   */
  void operator=(const StandardMCMCOutput &another);
  //MoveOutput* duplicate() const;
  
  /**
   * @brief Performs the push back operation.
   *
   * @param particle_in The particle.
   */
  void push_back(const Particle &particle_in);
  
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
   * @brief Returns the iteration counter pointer.
   *
   * @return The result.
   */
  size_t* get_iteration_counter_pointer();
  
  /**
   * @brief Performs the increment counter operation.
   */
  void increment_counter();
  
  /**
   * @brief Performs the reset counter operation.
   */
  void reset_counter();
  
  /**
   * @brief Performs the mcmc adapt operation.
   */
  void mcmc_adapt();
  
  Particle move(RandomNumberGenerator &rng,
                const Particle &particle) const;
  
  Particle subsample_move(RandomNumberGenerator &rng,
                          const Particle &particle) const;
  
  /**
   * @brief Returns @c true if the termination criterion has been met.
   *
   * @return The result.
   */
  bool terminate() const;
  
  /**
   * @brief Returns the mcmc.
   *
   * @return The result.
   */
  virtual MCMC* get_mcmc()=0;
  /**
   * @brief Returns the mcmc.
   *
   * @return The result.
   */
  virtual const MCMC* get_mcmc() const=0;
  
  /**
   * @brief Returns the current algorithm parameters.
   *
   * @return The result.
   */
  Parameters get_current_algorithm_parameters() const;
  
protected:
  
  //virtual void specific_mcmc_adapt()=0;
  
  /** @brief The iteration counter. */
  size_t iteration_counter;
  
  // Stored here.
  /** @brief The termination. */
  MCMCTermination* termination;
  
  /**
   * @brief Copies the state of another MCMCTermination into this object.
   *
   * @param another The MCMCTermination instance to copy from.
   */
  void make_copy(const StandardMCMCOutput &another);
  
  /** @brief The output. */
  std::deque<Particle> output;
  
  /** @brief The algorithm parameters. */
  std::deque<Parameters> algorithm_parameters;
  
};
}

#endif
