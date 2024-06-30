#ifndef ProposalKernel_H
#define ProposalKernel_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "kernel.h"
#include "particle.h"
//#include "ensemble_member.h"
#include "distributions.h"
#include "ilike_header.h"

class MCMCAdaptor;
class SMCAdaptor;
class Transform;

class SMCOutput;

class CompositeProposalKernel;
class CompositeIndependentProposalKernel;
class ReinforceGradientEstimator;
class DirectGradientEstimatorOutput;

class ProposalKernel : public Kernel
{

public:

  ProposalKernel();
  virtual ~ProposalKernel();
  
  ProposalKernel(const Parameters &proposal_parameters_in);

  ProposalKernel(const ProposalKernel &another);

  void operator=(const ProposalKernel &another);
  virtual ProposalKernel* proposal_kernel_duplicate() const=0;

  Particle move(RandomNumberGenerator &rng,
                const Particle &particle) const;
  /*
  Particle move(RandomNumberGenerator &rng,
                Particle &particle,
                const Parameters &conditioned_on_parameters) const;
  */
  
  Particle subsample_move(RandomNumberGenerator &rng,
                          const Particle &particle) const;
  
  /*
  Particle subsample_move(RandomNumberGenerator &rng,
                          Particle &particle,
                          const Parameters &conditioned_on_parameters) const;
  */
  Particle subsample_move(RandomNumberGenerator &rng,
                          const std::string &variable,
                          const Particle &particle) const;
  
  /*
  Particle subsample_move(RandomNumberGenerator &rng,
                          const std::string &variable,
                          Particle &particle,
                          const Parameters &conditioned_on_parameters) const;
  */
  
  /*
  Particle move(RandomNumberGenerator &rng,
                      const Index* index,
                      Particle &particle) const;
  Particle move(RandomNumberGenerator &rng,
                      const Index* index,
                      Particle &particle,
                      const Parameters &conditioned_on_parameters) const;
  //Particle subsample_move(RandomNumberGenerator &rng,
  //                              Particle &particle) const;
  Particle subsample_move(RandomNumberGenerator &rng,
                                const Index* index,
                                Particle &particle,
                                const Parameters &conditioned_on_parameters) const;
  //Particle subsample_move(RandomNumberGenerator &rng,
  //                              const std::string &variable,
  //                              Particle &particle) const;
  Particle subsample_move(RandomNumberGenerator &rng,
                                const Index* index,
                                const std::string &variable,
                                Particle &particle,
                                const Parameters &conditioned_on_parameters) const;
  */
  
  void ensemble_adapt(EnsembleKalmanOutput* current_state);
  
  void smc_adapt(SMCOutput* current_state);
  
  void mcmc_adapt(const Particle &current_particle,
                  size_t iteration_counter);
  
  //void use_transform(const Particle &particle);
  
  double evaluate_kernel(const Particle &proposed_particle,
                         const Particle &old_particle) const;
  double subsample_evaluate_kernel(const Particle &proposed_particle,
                                   const Particle &old_particle) const;
  
  //double evaluate_kernel_for_mh(const Particle &proposed_particle) const;
  //double subsample_evaluate_kernel_for_mh(const Particle &proposed_particle) const;
  
  /*
  double subsample_evaluate_kernel(Particle &proposed_particle,
                                   Particle &old_particle,
                                   const Parameters &conditioned_on_parameters) const;
  */
  
  arma::mat gradient_of_log(const std::string &variable,
                            const Particle &proposed_particle,
                            const Particle &old_particle);
  
  /*
  arma::mat gradient_of_log(const std::string &variable,
                                    Particle &proposed_particle,
                                    Particle &old_particle,
                                    const Parameters &conditioned_on_parameters);
  */
  
  arma::mat subsample_gradient_of_log(const std::string &variable,
                                      const Particle &proposed_particle,
                                      const Particle &old_particle);
  
  /*
  arma::mat subsample_gradient_of_log(const std::string &variable,
                                              Particle &proposed_particle,
                                              Particle &old_particle,
                                              const Parameters &conditioned_on_parameters);
  */
  
  virtual Parameters simulate(RandomNumberGenerator &rng,
                              const Particle &particle) const=0;
  
  /*
  virtual Parameters simulate(RandomNumberGenerator &rng,
                              Particle &particle,
                              const Parameters &conditioned_on_parameters) const=0;
  */
  
  virtual Parameters subsample_simulate(RandomNumberGenerator &rng,
                                        const Particle &particle) const=0;
  
  /*
  virtual Parameters subsample_simulate(RandomNumberGenerator &rng,
                                        Particle &particle,
                                        const Parameters &conditioned_on_parameters) const=0;
  */
  
  virtual Parameters subsample_simulate(RandomNumberGenerator &rng,
                                        const std::string &variable,
                                        const Particle &particle) const=0;
  
  /*
  virtual Parameters subsample_simulate(RandomNumberGenerator &rng,
                                        const std::string &variable,
                                        Particle &particle,
                                        const Parameters &conditioned_on_parameters) const=0;
  */
  
  Transform* get_transform() const;
  
  virtual void set_proposal_parameters(Parameters* proposal_parameters_in)=0;
  
  virtual GradientEstimatorOutput* simulate_gradient_estimator_output() const=0;
  
  virtual std::vector<const ProposalKernel*> get_proposals() const=0;
  
  virtual void set_index(Index* index_in)=0;
  
  virtual void set_index_if_null(Index* index_in)=0;
  
  int get_instance_index() const;
  
  virtual bool can_be_evaluated() const=0;
  
  virtual void set_data(Data* data_in)=0;
  
// Mh has its own parameters.
  // Stochastic has some weights.
  // MH has sim prop and eval prop, take in params. Use current value in acceptance, Set current value if accepted.
  // Proposal needs to call simulate in all llhdoutputs

protected:
  
  friend CompositeProposalKernel;
  friend CompositeIndependentProposalKernel;
  friend ReinforceGradientEstimator;
  friend DirectGradientEstimatorOutput;
  
  virtual double specific_evaluate_kernel(const Particle &proposed_particle,
                                          const Particle &old_particle) const=0;
  
  /*
  virtual double specific_evaluate_kernel(Particle &proposed_particle,
                                          Particle &old_particle,
                                          const Parameters &conditioned_on_parameters) const=0;
  */
  
  virtual double specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                    const Particle &old_particle) const=0;
  
  /*
  virtual double specific_subsample_evaluate_kernel(Particle &proposed_particle,
                                                    Particle &old_particle,
                                                    const Parameters &conditioned_on_parameters) const=0;
  */
   
  virtual arma::mat specific_gradient_of_log(const std::string &variable,
                                             const Particle &proposed_particle,
                                             const Particle &old_particle)=0;
  
  /*
  virtual arma::mat specific_gradient_of_log(const std::string &variable,
                                             Particle &proposed_particle,
                                             Particle &old_particle,
                                             const Parameters &conditioned_on_parameters)=0;
  */
  
  virtual arma::mat specific_subsample_gradient_of_log(const std::string &variable,
                                                       const Particle &proposed_particle,
                                                       const Particle &old_particle)=0;
  
  /*
  virtual arma::mat specific_subsample_gradient_of_log(const std::string &variable,
                                                       Particle &proposed_particle,
                                                       Particle &old_particle,
                                                       const Parameters &conditioned_on_parameters)=0;
  */

  void make_copy(const ProposalKernel &another);
  
  //EvaluateLogMCMCProposalPtr proposal_evaluate;
  
  //SimulateMCMCProposalPtr proposal_simulate;
  
  // stored here
  MCMCAdaptor* mcmc_adaptor;
  
  // stored here
  SMCAdaptor* smc_adaptor;
  
  // stored here
  Transform* transform;
  
  static int instance_counter; // Static member variable to track instances
  int instance_index; // Instance-specific index
  
  //Parameters proposal_parameters;

};

#endif
