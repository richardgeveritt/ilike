#ifndef ProposalKernel_H
#define ProposalKernel_H

#include <RcppArmadillo.h>
#include <memory>
using namespace Rcpp;

#include <vector>

#include "kernel.h"
#include "particle.h"
#include "distributions.h"
#include "ilike_header.h"
#include "proposal_target.h"

namespace ilike
{

class MCMCAdaptor;
struct AdaptationStep;
class GaussianProposalInfo;
struct PrincipalComponentsInfo;
class SMCAdaptor;
class Transform;

class SMCOutput;

class CompositeProposalKernel;
class CompositeIndependentProposalKernel;
class ReinforceGradientEstimator;
class DirectGradientEstimatorOutput;

// log N(x; mu, scale * Sigma) from Sigma^{-1} and log|Sigma|: the expression
// the Gaussian kernels used, dmvnorm_using_precomp(x, mu, (1/scale)*inv,
// dim*log(scale)+logdet), without building the scaled d x d inverse when the
// scale is 1 (multiplying by 1.0 is exact, so the value is the same).
inline double dmvnorm_scaled_precomp(const arma::colvec &x,
                                     const arma::colvec &mu,
                                     const arma::mat &inv,
                                     double scale,
                                     double logdet)
{
  double dim = double(mu.n_rows);
  if (scale == 1.0)
    return dmvnorm_using_precomp(x, mu, inv, dim*log(scale)+logdet);
  else
    return dmvnorm_using_precomp(x, mu, (1.0/scale)*inv, dim*log(scale)+logdet);
}

class ProposalKernel : public Kernel
{
  
public:

  bool has_mcmc_adaptor() const { return this->mcmc_adaptor!=NULL; }
  // adapting with until = and not yet stopped (C212)
  bool mcmc_warming_up() const;
  
  ProposalKernel();
  virtual ~ProposalKernel();
  
  ProposalKernel(const Parameters &proposal_parameters_in);
  
  ProposalKernel(const ProposalKernel &another);
  
  void operator=(const ProposalKernel &another);
  virtual ProposalKernel* proposal_kernel_duplicate() const=0;
  
  Particle move(RandomNumberGenerator &rng,
                const Particle &particle) const;

  // As above, but only the factors in new_outputs_index get new outputs in the
  // proposed particle; the others are copied from particle with their cached
  // values (see FactorVariables::copy_with_new_outputs). NULL means all
  // factors get new outputs (the behaviour of the overload above).
  // (virtual: the HMC proposal builds its particle with its momenta's
  // energies, C208)
  // (dev/plan.md 3.21) as a measurement model: the largest log density of
  // its variables given the others, over the others (+inf if not known)
  virtual double log_density_maximum(const Parameters &parameters) const { return arma::datum::inf; }

  virtual Particle move(RandomNumberGenerator &rng,
                        const Particle &particle,
                        const Index* new_outputs_index) const;

  // As move(rng, particle, new_outputs_index), but the proposal is built in
  // scratch, reusing its objects where possible (dev/plan.md, Stage 3 (a)).
  // Same random numbers and the same result as move.
  virtual void move_into(RandomNumberGenerator &rng,
                         const Particle &particle,
                         const Index* new_outputs_index,
                         Particle &scratch) const;
  
  Particle subsample_move(RandomNumberGenerator &rng,
                          const Particle &particle) const;
  
  Particle subsample_move(RandomNumberGenerator &rng,
                          const std::string &variable,
                          const Particle &particle) const;

  void ensemble_adapt(EnsembleKalmanOutput* current_state);
  
  void smc_adapt(SMCOutput* current_state);
  
  void mcmc_adapt(const Particle &current_particle,
                  size_t iteration_counter);

  // Adaptation (dev/plan.md 3.2): attaches an adaptor (owned) and points it
  // at this kernel; one adaptation step from a Metropolis(-Hastings) move.
  void set_mcmc_adaptor(MCMCAdaptor* adaptor_in);
  // adaptation from the SMC particles (dev/plan.md 3.2, B2): owned
  void set_smc_adaptor(SMCAdaptor* adaptor_in);
  SMCAdaptor* get_smc_adaptor() const { return this->smc_adaptor.get(); }
  void mcmc_adapt(const AdaptationStep &step);
  void prepare_mcmc_adaptation();
  void append_adaptation_state(std::vector<std::string> *names,
                               std::vector<double> *values,
                               bool final,
                               const std::string &prefix) const;

  // The proposal parameters of kernels written by the user (NULL otherwise).
  virtual Parameters* get_proposal_parameters() const;

  // Whether the Gaussian mean is part of the proposal (mirror), so that an
  // adaptor sets it to its running estimate of the target's mean.
  virtual bool proposal_uses_gaussian_mean() const;

  // The Gaussian parameters (mean, covariance, scale) of a variable, for
  // kernels that have them (random walk, Langevin, Barker); NULL otherwise.
  virtual GaussianProposalInfo* get_gaussian_proposal_info(const std::string &variable);

  // Per-component scales of a componentwise kernel (Algorithm 5); NULL otherwise.
  virtual arma::colvec* get_component_scales(const std::string &variable);

  // The principal directions of a principal-components kernel (Algorithm 8); NULL otherwise.
  virtual PrincipalComponentsInfo* get_principal_components(const std::string &variable);

  double evaluate_kernel(const Particle &proposed_particle,
                         const Particle &old_particle) const;
  // The same for two sets of parameters (e.g. a data model evaluated at the
  // data given a particle's state). Before this overload existed, such calls
  // converted both arguments to Particles; that is still what happens by
  // default (specific_evaluate_kernel_from_parameters), unless the kernel
  // evaluates the parameters directly.
  double evaluate_kernel(const Parameters &proposed_parameters,
                         const Parameters &old_parameters) const;
  double subsample_evaluate_kernel(const Particle &proposed_particle,
                                   const Particle &old_particle) const;

  arma::mat gradient_of_log(const std::string &variable,
                            const Particle &proposed_particle,
                            const Particle &old_particle);

  arma::mat subsample_gradient_of_log(const std::string &variable,
                                      const Particle &proposed_particle,
                                      const Particle &old_particle);

  virtual Parameters simulate(RandomNumberGenerator &rng,
                              const Particle &particle) const=0;

  // simulate, writing each proposed variable into target.matrix(name) instead
  // of returning a new Parameters (used by move_into when there is no
  // transform). Default: target.write_all(simulate(rng, particle)). Kernels
  // override it with the body of their simulate, writing into the target;
  // the random numbers and the values must be those of simulate.
  // (C269, PMMH P1) called at the start of a filter run with the run's
  // parameters, at a sequential point: a kernel may keep what depends only on
  // declared variables (see GaussianProposalFunctionsInfo::prepare); nothing
  // by default
  virtual void prepare_for_parameters(const Parameters &parameters) {}
  virtual void simulate_into(RandomNumberGenerator &rng,
                             const Particle &particle,
                             ProposalTarget &target) const;

  virtual Parameters subsample_simulate(RandomNumberGenerator &rng,
                                        const Particle &particle) const=0;

  virtual Parameters subsample_simulate(RandomNumberGenerator &rng,
                                        const std::string &variable,
                                        const Particle &particle) const=0;

  Transform* get_transform() const;
  // the scale the kernel proposes on (takes ownership; C209)
  void set_transform(Transform* transform_in);
  
  virtual void set_proposal_parameters(Parameters* proposal_parameters_in)=0;
  
  virtual GradientEstimatorOutput* simulate_gradient_estimator_output() const=0;
  
  virtual std::vector<const ProposalKernel*> get_proposals() const=0;
  // the kernels this proposal uses for a move from the particle (choose uses
  // one of its components, C224); by default all of them
  virtual std::vector<const ProposalKernel*> proposals_used(const Particle &particle) const { return this->get_proposals(); }
  
  virtual void set_index(Index* index_in)=0;
  
  virtual void set_index_if_null(Index* index_in)=0;
  
  int get_instance_index() const;
  // a particle's own tuning values for this kernel (esjd in SMC, C211), or
  // NULL: what they mean is the kernel's
  const double* tuning_values(const Particle &particle) const
  {
    return (particle.tuning != NULL) ? particle.tuning->for_kernel(this->get_instance_index()) : NULL;
  }
  
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
  
  // the steps of move and move_into after the simulation, for parameters
  // simulated by the kernel on its own space
  Particle proposed_from(const Particle &particle,
                         const Parameters &simulated,
                         const Index* new_outputs_index) const;
  void into_from(const Particle &particle,
                 const Parameters &simulated,
                 const Index* new_outputs_index,
                 Particle &scratch) const;
  void finish_into(const Particle &particle,
                   const Index* new_outputs_index,
                   Particle &scratch) const;

  virtual double specific_evaluate_kernel(const Particle &proposed_particle,
                                          const Particle &old_particle) const=0;

  // Called by evaluate_kernel(Parameters, Parameters) when the kernel has no
  // transform. Default: build the two Particles and call
  // specific_evaluate_kernel. Kernels whose density depends only on the
  // (untransformed) parameters override it to skip building the Particles.
  virtual double specific_evaluate_kernel_from_parameters(const Parameters &proposed_parameters,
                                                          const Parameters &old_parameters) const;

  virtual double specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                    const Particle &old_particle) const=0;

  virtual arma::mat specific_gradient_of_log(const std::string &variable,
                                             const Particle &proposed_particle,
                                             const Particle &old_particle)=0;

  virtual arma::mat specific_subsample_gradient_of_log(const std::string &variable,
                                                       const Particle &proposed_particle,
                                                       const Particle &old_particle)=0;

  void make_copy(const ProposalKernel &another);

  // stored here
  // owned (dev/plan.md 3.5, CH3)
  std::unique_ptr<MCMCAdaptor> mcmc_adaptor;
  
  // stored here
  // owned (dev/plan.md 3.5, CH3)
  std::unique_ptr<SMCAdaptor> smc_adaptor;
  
  // stored here
  // owned (dev/plan.md 3.5, CH3)
  std::unique_ptr<Transform> transform;
  
  static int instance_counter; // Static member variable to track instances
  int instance_index; // Instance-specific index

};
}

#endif
