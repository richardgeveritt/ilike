#ifndef PARAMETERPARTICLESIMULATOR_H
#define PARAMETERPARTICLESIMULATOR_H

#include "particle_simulator.h"

#include "function_pointers.h"
#include "simulation.h"

class ParameterParticleSimulator : public ParticleSimulator
{
public:

  ParameterParticleSimulator(SimulateDistributionPtr simulate_parameters_in);
  ParameterParticleSimulator(const ParameterParticleSimulator &another);
  virtual ~ParameterParticleSimulator(void);

  void operator=(const ParameterParticleSimulator &another);
  ParticleSimulator* duplicate() const;

  Particle operator()(RandomNumberGenerator &rng) const;

protected:

  void make_copy(const ParameterParticleSimulator &another);

  SimulateDistributionPtr simulate_parameters;
};

#endif
