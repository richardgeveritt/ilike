#ifndef PARTICLESIMULATOR_H
#define PARTICLESIMULATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <string>

#include "distributions.h"
#include "particle.h"

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
                            Factors* factors,
                            std::vector<ProposalKernel*>* proposals_to_transform_for,
                            std::vector<ProposalKernel*>* proposals_to_find_gradient_for) const=0;
  
  virtual Particle subsample_simulate(RandomNumberGenerator &rng,
                                      Factors* factors,
                                      std::vector<ProposalKernel*>* proposals_to_transform_for,
                                      std::vector<ProposalKernel*>* proposals_to_find_gradient_for) const=0;
  
  virtual Particle simulate(RandomNumberGenerator &rng,
                            Factors* factors,
                            std::vector<ProposalKernel*>* proposals_to_transform_for,
                            std::vector<ProposalKernel*>* proposals_to_find_gradient_for,
                            const Parameters &sequencer_parameters) const=0;
  
  virtual Particle subsample_simulate(RandomNumberGenerator &rng,
                                      Factors* factors,
                                      std::vector<ProposalKernel*>* proposals_to_transform_for,
                                      std::vector<ProposalKernel*>* proposals_to_find_gradient_for,
                                      const Parameters &sequencer_parameters) const=0;
  
  virtual Particle simulate(RandomNumberGenerator &rng,
                            Factors* factors,
                            std::vector<ProposalKernel*>* proposals_to_transform_for,
                            std::vector<ProposalKernel*>* proposals_to_find_gradient_for,
                            const Parameters &conditioned_on_parameters,
                            const Parameters &sequencer_parameters) const=0;
  
  virtual Particle subsample_simulate(RandomNumberGenerator &rng,
                                      Factors* factors,
                                      std::vector<ProposalKernel*>* proposals_to_transform_for,
                                      std::vector<ProposalKernel*>* proposals_to_find_gradient_for,
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
  virtual double subsample_evaluate(const Particle &input) const=0;
  
  /*
  virtual double evaluate(Particle &input,
                            const Parameters &conditioned_on_parameters) const=0;
  
  virtual double subsample_evaluate(Particle &input,
                                    const Parameters &conditioned_on_parameters) const=0;
  */

protected:
  
  std::string resample_variable_name;

  void make_copy(const ParticleSimulator &another);

};

#endif
