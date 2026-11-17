#ifndef SMCMARGINAL_H
#define SMCMARGINAL_H

#include "smc.h"
#include "proposal_kernel.h"

namespace ilike
{

class SMCOutput;
class MoveOutput;
class Index;

class SMCMarginal : public SMC
{
public:
  
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
  
  SMCMarginal(const SMCMarginal &another);
  virtual ~SMCMarginal();
  
  void operator=(const SMCMarginal &another);
  SMC* smc_duplicate() const;
  LikelihoodEstimator* duplicate() const;
  
  MoveOutput* move(RandomNumberGenerator &rng,
                   const Particle &particle) const;

  void weight_for_adapting_sequence(const Index* index,
                                    Particles &current_particles);

  MoveOutput* subsample_move(RandomNumberGenerator &rng,
                             const Particle &particle) const;

  void subsample_weight_for_adapting_sequence(const Index* index,
                                              Particles &current_particles);

protected:
  
  SMCOutput* specific_run();
  
  SMCOutput* specific_initialise_smc();
  void simulate_smc(SMCOutput* simulation);
  
  void evaluate_smc(SMCOutput* simulation);
  void evaluate_smcfixed_part_smc(SMCOutput* simulation);
  void evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation);
  
  void mcmc_move(SMCOutput* current_state);
  
  SMCOutput* specific_run(const Parameters &parameters);
  SMCOutput* specific_initialise_smc(const Parameters &parameters);
  void simulate_smc(SMCOutput* simulation,
                    const Parameters &conditioned_on_parameters);
  
  void evaluate_smc(SMCOutput* simulation,
                    const Parameters &conditioned_on_parameters);
  void evaluate_smcfixed_part_smc(SMCOutput* simulation,
                                  const Parameters &conditioned_on_parameters);
  void evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation,
                                                    const Parameters &conditioned_on_parameters);

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

  void make_copy(const SMCMarginal &another);
  
  // Stored here.
  std::unique_ptr<ProposalKernel> proposal_kernel;
  
  // stored here
  std::unique_ptr<Index> index;

};
}

#endif
