#ifndef PARTICLESIMULATOR_H
#define PARTICLESIMULATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <string>

#include "distributions.h"
#include "particle.h"

namespace ilike
{
  /**
   * @file particle_simulator.h
   * @brief Defines the Factors class.
   *
   * Provides factors functionality.
   *
   * @namespace ilike
   * @class Factors
   * @brief The factors class.
   */


class Factors;

class ParticleSimulator
{
public:
  
  /**
   * @brief Performs the particlesimulator operation.
   */
  ParticleSimulator();
  /**
   * @brief Performs the particlesimulator operation.
   *
   * @param resample_variable_name_in The resample variable name.
   */
  ParticleSimulator(const std::string &resample_variable_name_in);
  /**
   * @brief Performs the particlesimulator operation.
   *
   * @param another The Factors instance to copy from.
   */
  ParticleSimulator(const ParticleSimulator &another);
  /**
   * @brief Performs the ~particlesimulator operation.
   */
  virtual ~ParticleSimulator();
  
  /**
   * @brief Assignment operator for Factors.
   *
   * @param another The Factors instance to copy from.
   */
  void operator=(const ParticleSimulator &another);
  /**
   * @brief Creates a deep copy of this Factors object.
   *
   * @return The result.
   */
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
  
  /*
   virtual void simulate_and_transform(RandomNumberGenerator &rng,
   Particle* new_particle,
   Factors* factors,
   Transform* transform,
   bool store_raw) const=0;
   virtual void simulate_and_transform(RandomNumberGenerator &rng,
   Particle* new_particle,
   Factors* factors,
   Transform* transform,
   bool store_raw,
   const Parameters &conditioned_on_parameters) const=0;
   
   virtual void subsample_simulate_and_transform(RandomNumberGenerator &rng,
   Particle* new_particle,
   Factors* factors,
   Transform* transform,
   bool store_raw,
   const Parameters &conditioned_on_parameters) const=0;
   */
  
  virtual double evaluate(const Particle &input) const=0;
  /**
   * @brief Performs the subsample evaluate operation.
   *
   * @param input The input.
   *
   * @return The result.
   */
  virtual double subsample_evaluate(const Particle &input) const=0;
  
  /*
   virtual double evaluate(Particle &input,
   const Parameters &conditioned_on_parameters) const=0;
   
   virtual double subsample_evaluate(Particle &input,
   const Parameters &conditioned_on_parameters) const=0;
   */
  
protected:
  
  /** @brief The resample variable name. */
  std::string resample_variable_name;
  
  /**
   * @brief Copies the state of another Factors into this object.
   *
   * @param another The Factors instance to copy from.
   */
  void make_copy(const ParticleSimulator &another);
  
};
}

#endif
