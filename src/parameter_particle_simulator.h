#ifndef PARAMETERPARTICLESIMULATOR_H
#define PARAMETERPARTICLESIMULATOR_H

#include "particle_simulator.h"

#include "function_pointers.h"
#include "distributions.h"

class LikelihoodEstimator;

class ParameterParticleSimulator : public ParticleSimulator
{
public:

  ParameterParticleSimulator();
  ParameterParticleSimulator(SimulateDistributionPtr simulate_parameters_in,
                             const std::vector<LikelihoodEstimator*> &likelihood_estimators_in);
  virtual ~ParameterParticleSimulator(void);

  ParameterParticleSimulator(const ParameterParticleSimulator &another);
  void operator=(const ParameterParticleSimulator &another);
  ParticleSimulator* duplicate() const;

  Particle operator()(RandomNumberGenerator &rng);

protected:

  void make_copy(const ParameterParticleSimulator &another);

  SimulateDistributionPtr simulate_parameters;

  std::vector<LikelihoodEstimator*> likelihood_estimators;
};

#endif
