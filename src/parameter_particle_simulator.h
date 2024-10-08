#ifndef PARAMETERPARTICLESIMULATOR_H
#define PARAMETERPARTICLESIMULATOR_H

#include "particle_simulator.h"

#include "ilike_header.h"
#include "distributions.h"

namespace ilike
{
class IndependentProposalKernel;
class LikelihoodEstimator;

class ParameterParticleSimulator : public ParticleSimulator
{
public:
  
  ParameterParticleSimulator();
  ParameterParticleSimulator(const IndependentProposalKernel* proposal_in,
                             const std::vector<LikelihoodEstimator*> &likelihood_estimators_in);
  //ParameterParticleSimulator(SimulateDistributionPtr simulate_parameters_in,
  //                           const std::vector<LikelihoodEstimator*> &likelihood_estimators_in,
  //                          const std::string &resample_variable_name_in);
  virtual ~ParameterParticleSimulator();
  
  ParameterParticleSimulator(const ParameterParticleSimulator &another);
  void operator=(const ParameterParticleSimulator &another);
  ParticleSimulator* duplicate() const;
  
  Particle simulate(RandomNumberGenerator &rng,
                    const Factors* factors,
                    const std::vector<const ProposalKernel*>* proposals_to_transform_for,
                    const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for) const;
  
  Particle subsample_simulate(RandomNumberGenerator &rng,
                              const Factors* factors,
                              const std::vector<const ProposalKernel*>* proposals_to_transform_for,
                              const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for) const;
  
  Particle simulate(RandomNumberGenerator &rng,
                    const Factors* factors,
                    const std::vector<const ProposalKernel*>* proposals_to_transform_for,
                    const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for,
                    const Parameters &sequencer_parameters) const;
  
  Particle subsample_simulate(RandomNumberGenerator &rng,
                              const Factors* factors,
                              const std::vector<const ProposalKernel*>* proposals_to_transform_for,
                              const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for,
                              const Parameters &sequencer_parameters) const;
  
  Particle simulate(RandomNumberGenerator &rng,
                    const Factors* factors,
                    const std::vector<const ProposalKernel*>* proposals_to_transform_for,
                    const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for,
                    const Parameters &conditioned_on_parameters,
                    const Parameters &sequencer_parameters) const;
  
  Particle subsample_simulate(RandomNumberGenerator &rng,
                              const Factors* factors,
                              const std::vector<const ProposalKernel*>* proposals_to_transform_for,
                              const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for,
                              const Parameters &conditioned_on_parameters,
                              const Parameters &sequencer_parameters) const;
  
  
  /*
   void simulate_and_transform(RandomNumberGenerator &rng,
   Particle* new_particle,
   Factors* factors,
   Transform* transform,
   bool store_raw) const;
   void simulate_and_transform(RandomNumberGenerator &rng,
   Particle* new_particle,
   Factors* factors,
   Transform* transform,
   bool store_raw,
   const Parameters &conditioned_on_parameters) const;
   
   void subsample_simulate_and_transform(RandomNumberGenerator &rng,
   Particle* new_particle,
   Factors* factors,
   Transform* transform,
   bool store_raw,
   const Parameters &conditioned_on_parameters) const;
   */
  
  double evaluate(const Particle &input) const;
  
  double subsample_evaluate(const Particle &input) const;
  
  /*
   double evaluate(Particle &input,
   const Parameters &conditioned_on_parameters) const;
   
   double subsample_evaluate(Particle &input,
   const Parameters &conditioned_on_parameters) const;
   */
  
protected:
  
  void make_copy(const ParameterParticleSimulator &another);
  
  //SimulateDistributionPtr simulate_parameters;
  
  // stored here
  const IndependentProposalKernel* proposal;
  
  // Not stored here.
  std::vector<LikelihoodEstimator*> likelihood_estimators;
};
}

#endif
