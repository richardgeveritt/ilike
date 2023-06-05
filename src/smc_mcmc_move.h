#ifndef SMCMCMCMOVE_H
#define SMCMCMCMOVE_H

#include "smc.h"
#include "ilike_header.h"
#include "independent_proposal_kernel.h"

class SMCOutput;
class MoveOutput;
class Index;

class SMCMCMCMove : public SMC
{
public:

  SMCMCMCMove();
  
  // Multiple MCMC chains.
  SMCMCMCMove(RandomNumberGenerator* rng_in,
              size_t* seed_in,
              Data* data_in,
              const Parameters &algorithm_parameters,
              size_t lag_in,
              size_t lag_proposed_in,
              MCMC* mcmc_in,
              EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
              EvaluateLogDistributionPtr evaluate_log_prior_in,
              const std::vector<Parameters> &initial_points_in,
              const arma::colvec &log_probabilities_of_initial_values_in,
              bool parallel_in,
              size_t grain_size_in,
              const std::string &results_name_in);
  
  // Multiple MCMC chains.
  SMCMCMCMove(RandomNumberGenerator* rng_in,
              size_t* seed_in,
              Data* data_in,
              const Parameters &algorithm_parameters,
              size_t lag_in,
              size_t lag_proposed_in,
              MCMC* mcmc_in,
              const std::vector<LikelihoodEstimator*> &likelihood_estimators_in,
              const std::vector<Parameters> &initial_points_in,
              const arma::colvec &log_probabilities_of_initial_values_in,
              bool parallel_in,
              size_t grain_size_in,
              const std::string &results_name_in);
  
  SMCMCMCMove(RandomNumberGenerator* rng_in,
              size_t* seed_in,
              Data* data_in,
              size_t number_of_particles_in,
              size_t lag_in,
              size_t lag_proposed_in,
              MCMC* mcmc_in,
              EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
              EvaluateLogDistributionPtr evaluate_log_prior_in,
              SimulateDistributionPtr simulate_proposal_in,
              EvaluateLogDistributionPtr evaluate_log_proposal_in,
              bool mcmc_at_last_step_in,
              bool parallel_in,
              size_t grain_size_in,
              const std::string &results_name_in);

  // SMC with annealing.
  SMCMCMCMove(RandomNumberGenerator* rng_in,
              size_t* seed_in,
              Data* data_in,
              size_t number_of_particles_in,
              size_t lag_in,
              size_t lag_proposed_in,
              MCMC* mcmc_in,
              double resampling_desired_ess_in,
              double annealing_desired_cess_in,
              size_t number_of_bisections_in,
              EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
              EvaluateLogDistributionPtr evaluate_log_prior_in,
              SimulateDistributionPtr simulate_proposal_in,
              EvaluateLogDistributionPtr evaluate_log_proposal_in,
              bool mcmc_at_last_step_in,
              bool parallel_in,
              size_t grain_size_in,
              const std::string &results_name_in);
  
  // SMC with annealing.
  SMCMCMCMove(RandomNumberGenerator* rng_in,
              size_t* seed_in,
              Data* data_in,
              size_t number_of_particles_in,
              size_t lag_in,
              size_t lag_proposed_in,
              MCMC* mcmc_in,
              double resampling_desired_ess_in,
              double annealing_desired_cess_in,
              size_t number_of_bisections_in,
              const std::vector<double> &temperatures_in,
              EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
              EvaluateLogDistributionPtr evaluate_log_prior_in,
              SimulateDistributionPtr simulate_proposal_in,
              EvaluateLogDistributionPtr evaluate_log_proposal_in,
              bool mcmc_at_last_step_in,
              bool parallel_in,
              size_t grain_size_in,
              const std::string &results_name_in);
  
  SMCMCMCMove(RandomNumberGenerator* rng_in,
              size_t* seed_in,
              Data* data_in,
              size_t number_of_particles_in,
              size_t lag_in,
              size_t lag_proposed_in,
              MCMC* mcmc_in,
              double resampling_desired_ess_in,
              double annealing_desired_cess_in,
              size_t number_of_bisections_in,
              SMCTermination* termination_in,
              const std::string &sequence_variable_in,
              const std::vector<double> &schedule_in,
              const std::vector<LikelihoodEstimator*> &likelihood_estimators_in,
              IndependentProposalKernel* proposal_in,
              bool proposal_is_evaluated_in,
              bool smcfixed_flag_in,
              bool sequencer_limit_is_fixed_in,
              bool mcmc_at_last_step_in,
              bool parallel_in,
              size_t grain_size_in,
              const std::string &results_name_in);
  
  SMCMCMCMove(const SMCMCMCMove &another);
  
  virtual ~SMCMCMCMove(void);

  void operator=(const SMCMCMCMove &another);
  SMC* smc_duplicate() const;
  LikelihoodEstimator* duplicate() const;
  
  MoveOutput* move(RandomNumberGenerator &rng,
                   Particle &particle);
  
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
                             Particle &particle);
  
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
  
  /*
  void set_multiple_mcmc(RandomNumberGenerator* rng_in,
                         size_t* seed_in,
                         const Data* data_in,
                         MCMC* mcmc_in,
                         SimulateDistributionPtr simulate_proposal_in,
                         EvaluateLogDistributionPtr evaluate_log_prior_in,
                         EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
                         bool parallel_in,
                         size_t grain_size_in);
  */
  
  SMCOutput* specific_run();
  
  //SMCOutput* specific_run(const std::string &directory_name);
  
  SMCOutput* specific_initialise_smc();
  void simulate_smc(SMCOutput* simulation);
  void evaluate_smc(SMCOutput* simulation);
  void evaluate_smcfixed_part_smc(SMCOutput* simulation);
  void evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation);
  
  void mcmc_move(SMCOutput* current_state);
  
  SMCOutput* specific_run(const Parameters &parameters);
  
  //SMCOutput* specific_run(const std::string &directory_name,
  //                        const Parameters &parameters);
  
  SMCOutput* specific_initialise_smc(const Parameters &parameters);
  void simulate_smc(SMCOutput* simulation,
                    const Parameters &conditioned_on_parameters);
  
  void evaluate_smc(SMCOutput* simulation,
                    const Parameters &conditioned_on_parameters);
  void evaluate_smcfixed_part_smc(SMCOutput* simulation,
                                  const Parameters &conditioned_on_parameters);
  void evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation,
                                                    const Parameters &conditioned_on_parameters);
  
  /*
  void mcmc_move(SMCOutput* current_state,
                 const Parameters &conditioned_on_parameters);
  */
  
  void subsample_simulate_smc(SMCOutput* simulation);
  
  void subsample_evaluate_smc(SMCOutput* simulation);
  void subsample_evaluate_smcfixed_part_smc(SMCOutput* simulation);
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

  void make_copy(const SMCMCMCMove &another);
  
  // Stored here.
  MCMC* mcmc;
  
  // stored here
  Index* index;
  
  bool mcmc_at_last_step;

  //void smc_step();

  //void weight_update();
};

#endif
