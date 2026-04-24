#ifndef SMCOUTPUT_H
#define SMCOUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <deque>
#include <chrono>

#include "likelihood_estimator_output.h"
#include "particles.h"
#include "sequencer.h"

namespace ilike
{
  /**
   * @file smc_output.h
   * @brief Defines the SMC class.
   *
   * Provides smc functionality.
   *
   * @namespace ilike
   * @class SMC
   * @brief The smc class.
   */


class SMC;
class ImportanceSampler;
class SMCMCMCMove;
class SMCMarginal;
class SMCGeneric;
class ParticleFilter;

class SMCOutput : public LikelihoodEstimatorOutput
{
  
public:
  
  /**
   * @brief Performs the smcoutput operation.
   */
  SMCOutput();
  /**
   * @brief Performs the ~smcoutput operation.
   */
  virtual ~SMCOutput();
  
  SMCOutput(SMC* estimator_in,
            size_t lag_in,
            size_t lag_proposed_in,
            const std::string &results_name_in);
  
  /**
   * @brief Performs the smcoutput operation.
   *
   * @param another The SMC instance to copy from.
   */
  SMCOutput(const SMCOutput &another);
  /**
   * @brief Assignment operator for SMC.
   *
   * @param another The SMC instance to copy from.
   */
  void operator=(const SMCOutput &another);
  /**
   * @brief Creates a deep copy of this SMC object.
   *
   * @return The result.
   */
  LikelihoodEstimatorOutput* duplicate() const;
  /**
   * @brief Creates a deep copy and returns it as a smc pointer.
   *
   * @return The result.
   */
  SMCOutput* smc_duplicate() const;
  
  // Additional functions compared to some other LikelihoodEstimators, since we often run SMC without it simply being to estimate a likelihood.
  /**
   * @brief Simulates the required variables.
   */
  void simulate();
  /**
   * @brief Evaluates the smcfixed part.
   */
  void evaluate_smcfixed_part();
  /**
   * @brief Evaluates the smcadaptive part given smcfixed.
   */
  void evaluate_smcadaptive_part_given_smcfixed();
  
  /**
   * @brief Simulates the required variables.
   *
   * @param parameters The parameters.
   */
  void simulate(const Parameters &parameters);
  /**
   * @brief Evaluates the smcfixed part.
   *
   * @param parameters The parameters.
   */
  void evaluate_smcfixed_part(const Parameters &parameters);
  /**
   * @brief Evaluates the smcadaptive part given smcfixed.
   *
   * @param parameters The parameters.
   */
  void evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters);
  
  /**
   * @brief Performs the subsample simulate operation.
   */
  void subsample_simulate();
  /**
   * @brief Performs the subsample evaluate smcfixed part operation.
   */
  void subsample_evaluate_smcfixed_part();
  /**
   * @brief Performs the subsample evaluate smcadaptive part given smcfixed operation.
   */
  void subsample_evaluate_smcadaptive_part_given_smcfixed();
  
  /**
   * @brief Performs the subsample simulate operation.
   *
   * @param parameters The parameters.
   */
  void subsample_simulate(const Parameters &parameters);
  /**
   * @brief Performs the subsample evaluate smcfixed part operation.
   *
   * @param parameters The parameters.
   */
  void subsample_evaluate_smcfixed_part(const Parameters &parameters);
  /**
   * @brief Performs the subsample evaluate smcadaptive part given smcfixed operation.
   *
   * @param parameters The parameters.
   */
  void subsample_evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters);
  
  /**
   * @brief Performs the back operation.
   *
   * @return The result.
   */
  Particles back() const;
  /**
   * @brief Performs the back operation.
   *
   * @return The result.
   */
  Particles& back();
  
  /**
   * @brief Performs the end operation.
   *
   * @return The result.
   */
  std::deque<Particles>::iterator end();
  /**
   * @brief Performs the end operation.
   *
   * @return The result.
   */
  std::deque<Particles>::const_iterator end() const;
  
  /**
   * @brief Performs the calculate latest log normalising constant ratio operation.
   *
   * @return The result.
   */
  double calculate_latest_log_normalising_constant_ratio();
  
  arma::mat get_gradient_of_log(const std::string &variable,
                                const Parameters &x);
  arma::mat subsample_get_gradient_of_log(const std::string &variable,
                                          const Parameters &x);
  
  /**
   * @brief Performs the add particles operation.
   *
   * @return The result.
   */
  Particles* add_particles();
  /**
   * @brief Performs the add particles operation.
   *
   * @param most_recent_particles The most recent particles.
   *
   * @return The result.
   */
  Particles* add_particles(Particles* most_recent_particles);
  /**
   * @brief Performs the add proposed particles operation.
   *
   * @param latest_proposals The latest proposals.
   */
  void add_proposed_particles(const Particles &latest_proposals);
  //void initialise_unnormalised_log_incremental_weights(const arma::colvec &latest_unnormalised_log_incremental_weights);
  //void initialise_next_step();
  /**
   * @brief Updates weights.
   *
   * @param latest_unnormalised_log_incremental_weights The latest unnormalised log incremental weights.
   */
  void update_weights(const arma::colvec &latest_unnormalised_log_incremental_weights);
  /**
   * @brief Performs the normalise and resample weights operation.
   */
  void normalise_and_resample_weights();
  /**
   * @brief Resamples the particle population.
   */
  void resample();
  /**
   * @brief Performs the mcmc move operation.
   */
  void mcmc_move();
  
  /**
   * @brief Returns the likelihood estimator.
   *
   * @return The result.
   */
  LikelihoodEstimator* get_likelihood_estimator() const;
  
  /**
   * @brief Performs the number of smc iterations operation.
   *
   * @return The result.
   */
  size_t number_of_smc_iterations() const;
  
  /**
   * @brief Performs the increment smc iteration operation.
   */
  void increment_smc_iteration();
  /**
   * @brief Performs the decrement smc iteration operation.
   */
  void decrement_smc_iteration();
  
  /**
   * @brief Sets the time.
   */
  void set_time();
  /**
   * @brief Sets the llhd.
   *
   * @param llhd_in The llhd.
   */
  void set_llhd(double llhd_in);
  
  /**
   * @brief Sets the time and reset start.
   */
  void set_time_and_reset_start();
  
  /**
   * @brief Prints the object's state to an output stream.
   *
   * @param os The os.
   */
  void print(std::ostream &os) const;
  
  double log_likelihood_pre_last_step;
  
  /**
   * @brief Closes any open file streams.
   */
  void close_ofstreams();
  
  /**
   * @brief Performs the forget you were already written to file operation.
   */
  void forget_you_were_already_written_to_file();
  
  /**
   * @brief Returns @c true if the termination criterion has been met.
   */
  void terminate();
  
  std::deque<double> times;
  
  std::deque<double> llhds;
  
protected:
  
  //friend ImportanceSampler;
  //friend SMCMCMCMove;
  //friend SMC;
  //friend SMCMarginal;
  //friend SMCGeneric;
  //friend Sequencer;
  //friend ParticleFilter;
  
  // Stored in ModelAndAlgorithm or in main.
  /** @brief The estimator. */
  SMC* estimator;
  
  /** @brief The all particles. */
  std::deque<Particles> all_particles;
  
  /** @brief The all proposed. */
  std::deque<Particles> all_proposed;
  
  /** @brief The lag. */
  size_t lag;
  /** @brief The output lag. */
  size_t output_lag;
  
  /** @brief The lag proposed. */
  size_t lag_proposed;
  
  /** @brief The smc iteration. */
  size_t smc_iteration;
  
  /** @brief The iteration written to file. */
  int iteration_written_to_file;
  
  /** @brief The results name. */
  std::string results_name;
  
  /** @brief The start time. */
  std::chrono::high_resolution_clock::time_point start_time;
  
  void write_to_file(const std::string &directory_name,
                     const std::string &index = "");
  
  /**
   * @brief Closes any open file streams.
   *
   * @param deque_index The deque index.
   */
  void close_ofstreams(size_t deque_index);
  
  /**
   * @brief Copies the state of another SMC into this object.
   *
   * @param another The SMC instance to copy from.
   */
  void make_copy(const SMCOutput &another);
  
};
}

#endif
