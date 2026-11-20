#ifndef CHOICEPROPOSALKERNEL_H
#define CHOICEPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <memory>
#include <vector>

#include "proposal_kernel.h"
#include "esjd_adaptor.h"

namespace ilike
{

// choose([q_1, ..., q_K]) in a Metropolis-Hastings move in SMC (Fearnhead &
// Taylor 2013, the kernel index as part of the tuning h; dev/plan.md 3.19
// stage S4, C224). Each particle carries an index k (Particle::tuning, the
// ChoiceAdaptor's values) and moves with q_k for the whole move: proposals,
// and both directions of the Metropolis-Hastings ratio, are q_k's, so the
// move is q_k's Metropolis-Hastings kernel, which leaves the target
// unchanged whatever k is. The components keep their own parameters and
// adaptors (esjd, pretune, gaussian); only the kernel a particle used scores
// its adaptor (proposals_used).
class ChoiceProposalKernel : public ProposalKernel
{
public:

  ChoiceProposalKernel();
  ChoiceProposalKernel(const std::vector<ProposalKernel*> &components_in); // (takes ownership)
  virtual ~ChoiceProposalKernel();

  ChoiceProposalKernel(const ChoiceProposalKernel &another);
  void operator=(const ChoiceProposalKernel &another);
  Kernel* duplicate() const;
  ProposalKernel* proposal_kernel_duplicate() const;

  Particle move(RandomNumberGenerator &rng,
                const Particle &particle,
                const Index* new_outputs_index) const;
  void move_into(RandomNumberGenerator &rng,
                 const Particle &particle,
                 const Index* new_outputs_index,
                 Particle &scratch) const;
  Parameters simulate(RandomNumberGenerator &rng,
                      const Particle &particle) const;
  Parameters subsample_simulate(RandomNumberGenerator &rng,
                                const Particle &particle) const;
  Parameters subsample_simulate(RandomNumberGenerator &rng,
                                const std::string &variable,
                                const Particle &particle) const;

  void smc_adapt(SMCOutput* current_state);
  void set_proposal_parameters(Parameters* proposal_parameters_in);
  GradientEstimatorOutput* simulate_gradient_estimator_output() const;
  std::vector<const ProposalKernel*> get_proposals() const;
  std::vector<const ProposalKernel*> proposals_used(const Particle &particle) const;
  void set_index(Index* index_in);
  void set_index_if_null(Index* index_in);
  bool can_be_evaluated() const;
  void set_data(Data* data_in);

  size_t number_of_components() const { return this->components.size(); }
  ProposalKernel* component(size_t k) const { return this->components[k].get(); }

protected:

  // the particle's component (its tuning value), checked
  const ProposalKernel* chosen(const Particle &particle) const;

  double specific_evaluate_kernel(const Particle &proposed_particle,
                                  const Particle &old_particle) const;
  double specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                            const Particle &old_particle) const;
  arma::mat specific_gradient_of_log(const std::string &variable,
                                     const Particle &proposed_particle,
                                     const Particle &old_particle);
  arma::mat specific_subsample_gradient_of_log(const std::string &variable,
                                               const Particle &proposed_particle,
                                               const Particle &old_particle);

  void make_copy(const ChoiceProposalKernel &another);
  std::vector<std::unique_ptr<ProposalKernel>> components;
};

// The kernel index of choose (C224): as esjd, a value per particle, scored
// by the mean of alpha |x~ - x|^2_M over the particle's Metropolis-Hastings
// iterations divided by the cost of the kernel it used (a random walk 1,
// langevin or barker 2, hmc its steps + 1, or cost = [...]), resampled with
// those weights; the perturbation moves an index to another kernel, chosen
// uniformly, with probability `change`. The first indices are uniform.
class ChoiceAdaptor : public ESJDAdaptor
{
public:

  ChoiceAdaptor();
  ChoiceAdaptor(const std::string &variable_in, const std::vector<double> &costs_in, double change_in, bool full_metric_in);
  virtual ~ChoiceAdaptor();

  ChoiceAdaptor(const ChoiceAdaptor &another);
  void operator=(const ChoiceAdaptor &another);
  SMCAdaptor* duplicate() const;

  bool is_choice() const { return true; }
  void score(size_t i, const Parameters &before, const Parameters &after, const ParticleTuning &tuning);
  void summary(std::vector<std::string> &names, std::vector<double> &numbers) const;

protected:

  double perturbed(RandomNumberGenerator &rng, size_t k, double v) const;
  void make_copy(const ChoiceAdaptor &another);
  std::vector<double> costs;
  double change;
};

}

#endif
