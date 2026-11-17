#ifndef INDEPENDENTPROPOSALKERNEL_H
#define INDEPENDENTPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "proposal_kernel.h"
#include "particle.h"
#include "distributions.h"
#include "ilike_header.h"

namespace ilike
{

class BarkerDynamicsProposalKernel;

class IndependentProposalKernel : public ProposalKernel
{
  
public:
  
  IndependentProposalKernel();
  virtual ~IndependentProposalKernel();
  
  IndependentProposalKernel(const IndependentProposalKernel &another);
  
  void operator=(const IndependentProposalKernel &another);
  virtual IndependentProposalKernel* independent_proposal_kernel_duplicate() const=0;
  
  // Mh has its own parameters.
  // Stochastic has some weights.
  // MH has sim prop and eval prop, take in params. Use current value in acceptance, Set current value if accepted.
  // Proposal needs to call simulate in all llhdoutputs
  
  Parameters simulate(RandomNumberGenerator &rng,
                      const Particle &particle) const;

  Parameters subsample_simulate(RandomNumberGenerator &rng,
                                const Particle &particle) const;

  Parameters subsample_simulate(RandomNumberGenerator &rng,
                                const std::string &variable,
                                const Particle &particle) const;

  virtual Parameters independent_simulate(RandomNumberGenerator &rng) const=0;
  virtual Parameters subsample_independent_simulate(RandomNumberGenerator &rng) const=0;
  virtual Parameters subsample_independent_simulate(RandomNumberGenerator &rng,
                                                    const std::string &variable) const=0;
  
  virtual Parameters independent_simulate(RandomNumberGenerator &rng,
                                          const Parameters &conditioned_on_parameters) const=0;

  // Append number_of_points draws to output, all with the same
  // conditioned_on_parameters (e.g. the simulations of a synthetic
  // likelihood at one parameter value). Same draws as calling
  // independent_simulate that many times; kernels override it to do the work
  // that does not depend on the draw only once.
  virtual void independent_simulate_many(RandomNumberGenerator &rng,
                                         const Parameters &conditioned_on_parameters,
                                         size_t number_of_points,
                                         std::vector<Parameters> &output) const
  {
    for (size_t i = 0; i < number_of_points; ++i)
      output.push_back(this->independent_simulate(rng, conditioned_on_parameters));
  }

  virtual Parameters subsample_independent_simulate(RandomNumberGenerator &rng,
                                                    const Parameters &conditioned_on_parameters) const=0;

  virtual Parameters subsample_independent_simulate(RandomNumberGenerator &rng,
                                                    const std::string &variable,
                                                    const Parameters &conditioned_on_parameters) const=0;

  virtual double evaluate_independent_kernel(const Parameters &proposed_particle) const=0;
  virtual double subsample_evaluate_independent_kernel(const Parameters &proposed_particle) const=0;
  
  virtual arma::mat independent_gradient_of_log(const std::string &variable,
                                                const Parameters &proposed_particle)=0;
  virtual arma::mat subsample_independent_gradient_of_log(const std::string &variable,
                                                          const Parameters &proposed_particle)=0;
  
protected:
  
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

  friend BarkerDynamicsProposalKernel;
  
  void make_copy(const IndependentProposalKernel &another);

};

}

#endif
