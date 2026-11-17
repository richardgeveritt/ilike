#ifndef PARTICLESIMULATOR_H
#define PARTICLESIMULATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <string>

#include "distributions.h"
#include "particle.h"

namespace ilike
{

class Factors;

class ParticleSimulator
{
public:
  
  ParticleSimulator();
  ParticleSimulator(const std::string &resample_variable_name_in);
  ParticleSimulator(const ParticleSimulator &another);
  virtual ~ParticleSimulator();
  
  void operator=(const ParticleSimulator &another);
  virtual ParticleSimulator* duplicate() const=0;
  
  virtual Particle simulate(RandomNumberGenerator &rng,
                            const Factors* factors,
                            const std::vector<const ProposalKernel*>* proposals_to_transform_for,
                            const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for) const=0;
  
  virtual Particle subsample_simulate(RandomNumberGenerator &rng,
                                      const Factors* factors,
                                      const std::vector<const ProposalKernel*>* proposals_to_transform_for,
                                      const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for) const=0;
  
  virtual Particle simulate(RandomNumberGenerator &rng,
                            const Factors* factors,
                            const std::vector<const ProposalKernel*>* proposals_to_transform_for,
                            const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for,
                            const Parameters &sequencer_parameters) const=0;
  
  virtual Particle subsample_simulate(RandomNumberGenerator &rng,
                                      const Factors* factors,
                                      const std::vector<const ProposalKernel*>* proposals_to_transform_for,
                                      const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for,
                                      const Parameters &sequencer_parameters) const=0;
  
  virtual Particle simulate(RandomNumberGenerator &rng,
                            const Factors* factors,
                            const std::vector<const ProposalKernel*>* proposals_to_transform_for,
                            const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for,
                            const Parameters &conditioned_on_parameters,
                            const Parameters &sequencer_parameters) const=0;
  
  virtual Particle subsample_simulate(RandomNumberGenerator &rng,
                                      const Factors* factors,
                                      const std::vector<const ProposalKernel*>* proposals_to_transform_for,
                                      const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for,
                                      const Parameters &conditioned_on_parameters,
                                      const Parameters &sequencer_parameters) const=0;

  // Simulate into an existing particle object: the result is that of
  // target = simulate(rng, factors, ..., sequencer_parameters), with the same
  // random numbers (dev/plan.md 2.6a). Default: exactly that.
  virtual void simulate_into(RandomNumberGenerator &rng,
                             Particle &target,
                             const Factors* factors,
                             const std::vector<const ProposalKernel*>* proposals_to_transform_for,
                             const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for,
                             const Parameters &sequencer_parameters) const
  {
    target = this->simulate(rng, factors, proposals_to_transform_for, proposals_to_find_gradient_for, sequencer_parameters);
  }

  // The same for simulate(rng, factors, ..., conditioned_on_parameters,
  // sequencer_parameters).
  virtual void simulate_into(RandomNumberGenerator &rng,
                             Particle &target,
                             const Factors* factors,
                             const std::vector<const ProposalKernel*>* proposals_to_transform_for,
                             const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for,
                             const Parameters &conditioned_on_parameters,
                             const Parameters &sequencer_parameters) const
  {
    target = this->simulate(rng, factors, proposals_to_transform_for, proposals_to_find_gradient_for,
                            conditioned_on_parameters, sequencer_parameters);
  }

  virtual double evaluate(const Particle &input) const=0;
  virtual double subsample_evaluate(const Particle &input) const=0;

protected:
  
  std::string resample_variable_name;
  
  void make_copy(const ParticleSimulator &another);
  
};
}

#endif
