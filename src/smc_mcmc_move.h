#ifndef SMCMCMCMOVE_H
#define SMCMCMCMOVE_H

#include "smc.h"
#include "function_pointers.h"

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
              size_t lag_in,
              size_t lag_proposed_in,
              MCMC* mcmc_in,
              const std::vector<Parameters> &initial_points_in,
              const arma::colvec &log_probabilities_of_initial_values_in,
              EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
              EvaluateLogDistributionPtr evaluate_log_prior_in,
              bool parallel_in,
              size_t grain_size_in);
  
  SMCMCMCMove(RandomNumberGenerator* rng_in,
              size_t* seed_in,
              Data* data_in,
              size_t number_of_particles_in,
              size_t lag_in,
              size_t lag_proposed_in,
              MCMC* mcmc_in,
              EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
              EvaluateLogDistributionPtr evaluate_log_prior_in,
              SimulateIndependentProposalPtr simulate_proposal_in,
              EvaluateLogDistributionPtr evaluate_log_proposal_in,
              bool parallel_in,
              size_t grain_size_in);

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
              EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
              EvaluateLogDistributionPtr evaluate_log_prior_in,
              SimulateIndependentProposalPtr simulate_proposal_in,
              EvaluateLogDistributionPtr evaluate_log_proposal_in,
              bool parallel_in,
              size_t grain_size_in);
  
  SMCMCMCMove(const SMCMCMCMove &another);
  
  virtual ~SMCMCMCMove(void);

  void operator=(const SMCMCMCMove &another);
  SMC* smc_duplicate() const;
  LikelihoodEstimator* duplicate() const;
  
  MoveOutput* move(RandomNumberGenerator &rng,
                            Particle &particle);
  
  void weight_for_adapting_sequence(Particles &current_particles);
  
  MoveOutput* move(RandomNumberGenerator &rng,
                            Particle &particle,
                            const Parameters &conditioned_on_parameters);
  
  void weight_for_adapting_sequence(Particles &current_particles,
                                    const Parameters &conditioned_on_parameters);
  
  MoveOutput* subsample_move(RandomNumberGenerator &rng,
                             Particle &particle,
                             const Parameters &conditioned_on_parameters);
  
  void subsample_weight_for_adapting_sequence(Particles &current_particles,
                                              const Parameters &conditioned_on_parameters);

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
  
  SMCOutput* initialise_smc();
  void simulate_smc(SMCOutput* simulation);
  void evaluate_smc(SMCOutput* simulation);
  void evaluate_smcfixed_part_smc(SMCOutput* simulation);
  void evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation);
  
  void mcmc_move(SMCOutput* current_state);
  
  
  SMCOutput* specific_run(const Parameters &parameters);
  
  SMCOutput* initialise_smc(const Parameters &parameters);
  void simulate_smc(SMCOutput* simulation,
                    const Parameters &conditioned_on_parameters);
  
  void evaluate_smc(SMCOutput* simulation,
                    const Parameters &conditioned_on_parameters);
  void evaluate_smcfixed_part_smc(SMCOutput* simulation,
                    const Parameters &conditioned_on_parameters);
  void evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation,
                                          const Parameters &conditioned_on_parameters);
  
  void mcmc_move(SMCOutput* current_state,
                 const Parameters &conditioned_on_parameters);
  
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
  
  // Stored here.
  SMCCriterion* smc_criterion;
  
  // stored here
  Index* index;

  //void smc_step();

  //void weight_update();
};

#endif
