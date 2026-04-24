#ifndef SMCMARGINAL_H
#define SMCMARGINAL_H

#include "smc.h"
#include "proposal_kernel.h"

namespace ilike
{
  /**
   * @file smc_marginal.h
   * @brief Defines the SMCOutput class.
   *
   * Stores and manages the output produced by SMC. Holds results such as log-likelihood estimates, samples, or diagnostics returned after running the associated algorithm.
   *
   * @namespace ilike
   * @class SMCOutput
   * @brief The smc output class.
   */


class SMCOutput;
class MoveOutput;
class Index;

class SMCMarginal : public SMC
{
public:
  
  /**
   * @brief Performs the smcmarginal operation.
   */
  SMCMarginal();
  
  SMCMarginal(RandomNumberGenerator* rng_in,
              size_t* seed_in,
              Data* data_in,
              size_t number_of_particles_in,
              size_t lag_in,
              size_t lag_proposed_in,
              const std::vector<double> &temperatures_in,
              ProposalKernel* proposal_kernel_in,
              EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
              EvaluateLogDistributionPtr evaluate_log_prior_in,
              SimulateDistributionPtr simulate_proposal_in,
              EvaluateLogDistributionPtr evaluate_log_proposal_in,
              bool transform_proposed_particles,
              bool parallel_in,
              size_t grain_size_in,
              const std::string &results_name_in);
  
  /**
   * @brief Performs the smcmarginal operation.
   *
   * @param another The SMCOutput instance to copy from.
   */
  SMCMarginal(const SMCMarginal &another);
  /**
   * @brief Performs the ~smcmarginal operation.
   */
  virtual ~SMCMarginal();
  
  /**
   * @brief Assignment operator for SMCOutput.
   *
   * @param another The SMCOutput instance to copy from.
   */
  void operator=(const SMCMarginal &another);
  /**
   * @brief Creates a deep copy and returns it as a smc pointer.
   *
   * @return The result.
   */
  SMC* smc_duplicate() const;
  /**
   * @brief Creates a deep copy of this SMCOutput object.
   *
   * @return The result.
   */
  LikelihoodEstimator* duplicate() const;
  
  MoveOutput* move(RandomNumberGenerator &rng,
                   const Particle &particle) const;
  
  //void weight_for_adapting_sequence(Particles &current_particles);
  
  /*
   MoveOutput* move(RandomNumberGenerator &rng,
   Particle &particle,
   const Parameters &conditioned_on_parameters);
   */
  
  void weight_for_adapting_sequence(const Index* index,
                                    Particles &current_particles);
  
  /*
   void weight_for_adapting_sequence(const Index* index,
   Particles &current_particles,
   const Parameters &conditioned_on_parameters);
   */
  
  MoveOutput* subsample_move(RandomNumberGenerator &rng,
                             const Particle &particle) const;
  
  /*
   MoveOutput* subsample_move(RandomNumberGenerator &rng,
   Particle &particle,
   const Parameters &conditioned_on_parameters);
   */
  
  void subsample_weight_for_adapting_sequence(const Index* index,
                                              Particles &current_particles);
  
  /*
   void subsample_weight_for_adapting_sequence(const Index* index,
   Particles &current_particles,
   const Parameters &conditioned_on_parameters);
   */
  
protected:
  
  /**
   * @brief Class-specific implementation for run.
   *
   * @return The result.
   */
  SMCOutput* specific_run();
  //SMCOutput* specific_run(const std::string &directory_name);
  
  /**
   * @brief Class-specific implementation for initialise smc.
   *
   * @return The result.
   */
  SMCOutput* specific_initialise_smc();
  /**
   * @brief Simulates smc.
   *
   * @param simulation The simulation.
   */
  void simulate_smc(SMCOutput* simulation);
  
  /**
   * @brief Evaluates the smc.
   *
   * @param simulation The simulation.
   */
  void evaluate_smc(SMCOutput* simulation);
  /**
   * @brief Evaluates the smcfixed part smc.
   *
   * @param simulation The simulation.
   */
  void evaluate_smcfixed_part_smc(SMCOutput* simulation);
  /**
   * @brief Evaluates the smcadaptive part given smcfixed smc.
   *
   * @param simulation The simulation.
   */
  void evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation);
  
  /**
   * @brief Performs the mcmc move operation.
   *
   * @param current_state The current state.
   */
  void mcmc_move(SMCOutput* current_state);
  
  /**
   * @brief Class-specific implementation for run.
   *
   * @param parameters The parameters.
   *
   * @return The result.
   */
  SMCOutput* specific_run(const Parameters &parameters);
  /*
   SMCOutput* specific_run(const std::string &directory_name,
   const Parameters &parameters);
   */
  SMCOutput* specific_initialise_smc(const Parameters &parameters);
  void simulate_smc(SMCOutput* simulation,
                    const Parameters &conditioned_on_parameters);
  
  void evaluate_smc(SMCOutput* simulation,
                    const Parameters &conditioned_on_parameters);
  void evaluate_smcfixed_part_smc(SMCOutput* simulation,
                                  const Parameters &conditioned_on_parameters);
  void evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation,
                                                    const Parameters &conditioned_on_parameters);
  
  //void mcmc_move(SMCOutput* current_state,
  //               const Parameters &conditioned_on_parameters);
  
  /**
   * @brief Performs the subsample simulate smc operation.
   *
   * @param simulation The simulation.
   */
  void subsample_simulate_smc(SMCOutput* simulation);
  
  /**
   * @brief Performs the subsample evaluate smc operation.
   *
   * @param simulation The simulation.
   */
  void subsample_evaluate_smc(SMCOutput* simulation);
  /**
   * @brief Performs the subsample evaluate smcfixed part smc operation.
   *
   * @param simulation The simulation.
   */
  void subsample_evaluate_smcfixed_part_smc(SMCOutput* simulation);
  /**
   * @brief Performs the subsample evaluate smcadaptive part given smcfixed smc operation.
   *
   * @param simulation The simulation.
   */
  void subsample_evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation);
  
  void subsample_simulate_smc(SMCOutput* simulation,
                              const Parameters &conditioned_on_parameters);
  
  void subsample_evaluate_smc(SMCOutput* simulation,
                              const Parameters &conditioned_on_parameters);
  void subsample_evaluate_smcfixed_part_smc(SMCOutput* simulation,
                                            const Parameters &conditioned_on_parameters);
  void subsample_evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation,
                                                              const Parameters &conditioned_on_parameters);
  
  //void smc_update(SMCOutput* current_state);
  
  /**
   * @brief Copies the state of another SMCOutput into this object.
   *
   * @param another The SMCOutput instance to copy from.
   */
  void make_copy(const SMCMarginal &another);
  
  // Stored here.
  /** @brief The proposal kernel. */
  ProposalKernel* proposal_kernel;
  
  // stored here
  /** @brief The index. */
  Index* index;
  
  //void smc_step();
  
  //void weight_update();
};
}

#endif
