#ifndef PARAMETERPARTICLESIMULATOR_H
#define PARAMETERPARTICLESIMULATOR_H

#include "particle_simulator.h"

#include "ilike_header.h"
#include "distributions.h"

namespace ilike
{
  /**
   * @file parameter_particle_simulator.h
   * @brief Defines the IndependentProposalKernel class.
   *
   * An independent independent proposal kernel. Proposes new parameter values independently of the current state by sampling from a independent distribution.
   *
   * @namespace ilike
   * @class IndependentProposalKernel
   * @brief The independent proposal kernel class.
   */


class IndependentProposalKernel;
class LikelihoodEstimator;

class ParameterParticleSimulator : public ParticleSimulator
{
public:
  
  /**
   * @brief Performs the parameterparticlesimulator operation.
   */
  ParameterParticleSimulator();
  ParameterParticleSimulator(const IndependentProposalKernel* proposal_in,
                             const std::vector<LikelihoodEstimator*> &likelihood_estimators_in);
  //ParameterParticleSimulator(SimulateDistributionPtr simulate_parameters_in,
  //                           const std::vector<LikelihoodEstimator*> &likelihood_estimators_in,
  //                          const std::string &resample_variable_name_in);
  /**
   * @brief Performs the ~parameterparticlesimulator operation.
   */
  virtual ~ParameterParticleSimulator();
  
  /**
   * @brief Performs the parameterparticlesimulator operation.
   *
   * @param another The IndependentProposalKernel instance to copy from.
   */
  ParameterParticleSimulator(const ParameterParticleSimulator &another);
  /**
   * @brief Assignment operator for IndependentProposalKernel.
   *
   * @param another The IndependentProposalKernel instance to copy from.
   */
  void operator=(const ParameterParticleSimulator &another);
  /**
   * @brief Creates a deep copy of this IndependentProposalKernel object.
   *
   * @return The result.
   */
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
  
  /**
   * @brief Performs the subsample evaluate operation.
   *
   * @param input The input.
   *
   * @return The result.
   */
  double subsample_evaluate(const Particle &input) const;
  
  /*
   double evaluate(Particle &input,
   const Parameters &conditioned_on_parameters) const;
   
   double subsample_evaluate(Particle &input,
   const Parameters &conditioned_on_parameters) const;
   */
  
protected:
  
  /**
   * @brief Copies the state of another IndependentProposalKernel into this object.
   *
   * @param another The IndependentProposalKernel instance to copy from.
   */
  void make_copy(const ParameterParticleSimulator &another);
  
  //SimulateDistributionPtr simulate_parameters;
  
  // stored here
  /** @brief The proposal. */
  const IndependentProposalKernel* proposal;
  
  // Not stored here.
  /** @brief The likelihood estimators. */
  std::vector<LikelihoodEstimator*> likelihood_estimators;
};
}

#endif
