#include "mcmc.h"
#include "iterations_mcmc_termination.h"
#include "move_output.h"
#include "standard_mcmc_output.h"
#include "mcmc_adaptor.h"
#include "stochastic_scan_mcmc.h"
#include "deterministic_scan_mcmc.h"
#include "index.h"

MCMC::MCMC()
  :Kernel()
{
  this->iteration_counter = 0;
}

MCMC::MCMC(size_t number_of_iterations_in)
  :Kernel()
{
  this->iteration_counter = 0;
  this->termination = new IterationsMCMCTermination(number_of_iterations_in,
                                                    &this->iteration_counter);
}

MCMC::~MCMC()
{
  if (this->termination!=NULL)
    delete this->termination;
  
  //if (this->adaptor!=NULL)
  //  delete this->adaptor;
}

MCMC::MCMC(const MCMC &another)
  :Kernel(another)
{
  this->make_copy(another);
}

void MCMC::operator=(const MCMC &another)
{
  if(this == &another)
    return;
  
  if (this->termination!=NULL)
    delete this->termination;
  
  //if (this->adaptor!=NULL)
  //  delete this->adaptor;

  Kernel::operator=(another);
  this->make_copy(another);
}

void MCMC::make_copy(const MCMC &another)
{
  if (another.termination!=NULL)
    this->termination = another.termination->duplicate();
  
  //if (another.adaptor!=NULL)
  //  this->adaptor = another.adaptor->duplicate();
  
  this->iteration_counter = another.iteration_counter;
}

MoveOutput* MCMC::run(RandomNumberGenerator &rng,
                      Particle &particle)

{
  StandardMCMCOutput* mcmc_output = new StandardMCMCOutput();
  Particle current_particle = particle;
  while (!this->termination->terminate())
  {
    current_particle = this->move(rng,
                                  current_particle);
    this->iteration_counter = this->iteration_counter + 1;
    this->mcmc_adapt(current_particle,
                     this->iteration_counter);
    mcmc_output->push_back(current_particle);
  }
  return mcmc_output;
}

/*
MoveOutput* MCMC::run(RandomNumberGenerator &rng,
                      Particle &particle,
                      const Parameters &conditioned_on_parameters)

{
  StandardMCMCOutput* mcmc_output = new StandardMCMCOutput();
  Particle current_particle = particle;
  while (!this->termination->terminate())
  {
    current_particle = this->move(rng,
                                  current_particle,
                                  conditioned_on_parameters);
    this->iteration_counter = this->iteration_counter + 1;
    this->mcmc_adapt(current_particle,
                     this->iteration_counter);
    mcmc_output->push_back(current_particle);
  }
  return mcmc_output;
}
*/

MoveOutput* MCMC::subsample_run(RandomNumberGenerator &rng,
                                Particle &particle)

{
  StandardMCMCOutput* mcmc_output = new StandardMCMCOutput();
  Particle current_particle = particle;
  while (!this->termination->terminate())
  {
    current_particle = this->subsample_move(rng,
                                            current_particle);
    this->mcmc_adapt(current_particle,
                     this->iteration_counter);
    this->iteration_counter = this->iteration_counter + 1;
    mcmc_output->push_back(current_particle);
  }
  return mcmc_output;
}

/*
MoveOutput* MCMC::subsample_run(RandomNumberGenerator &rng,
                                Particle &particle,
                                const Parameters &conditioned_on_parameters)

{
  StandardMCMCOutput* mcmc_output = new StandardMCMCOutput();
  Particle current_particle = particle;
  while (!this->termination->terminate())
  {
    current_particle = this->subsample_move(rng,
                                            current_particle,
                                            conditioned_on_parameters);
    this->mcmc_adapt(current_particle,
                     this->iteration_counter);
    this->iteration_counter = this->iteration_counter + 1;
    mcmc_output->push_back(current_particle);
  }
  return mcmc_output;
}
*/

/*
EnsembleMember MCMC::run(RandomNumberGenerator &rng,
                         EnsembleMember &particle)

{
  StandardMCMCOutput* mcmc_output = new StandardMCMCOutput();
  EnsembleMember current_particle = particle;
  while (!this->termination->terminate())
  {
    current_particle = this->move(rng,
                                  current_particle);
    this->iteration_counter = this->iteration_counter + 1;
    this->mcmc_adapt(&current_particle,
                     this->iteration_counter);
  }
  return current_particle;
}

EnsembleMember MCMC::run(RandomNumberGenerator &rng,
                         EnsembleMember &particle,
                         const Parameters &conditioned_on_parameters)

{
  StandardMCMCOutput* mcmc_output = new StandardMCMCOutput();
  EnsembleMember current_particle = particle;
  while (!this->termination->terminate())
  {
    current_particle = this->move(rng,
                                  current_particle,
                                  conditioned_on_parameters);
    this->iteration_counter = this->iteration_counter + 1;
    this->mcmc_adapt(&current_particle,
                     this->iteration_counter);
  }
  return current_particle;
}

EnsembleMember MCMC::subsample_run(RandomNumberGenerator &rng,
                                   EnsembleMember &particle,
                                   const Parameters &conditioned_on_parameters)

{
  StandardMCMCOutput* mcmc_output = new StandardMCMCOutput();
  EnsembleMember current_particle = particle;
  while (!this->termination->terminate())
  {
    current_particle = this->subsample_move(rng,
                                            current_particle,
                                            conditioned_on_parameters);
    this->mcmc_adapt(&current_particle,
                     this->iteration_counter);
    this->iteration_counter = this->iteration_counter + 1;
  }
  return current_particle;
}
*/

void MCMC::mcmc_adapt(Particle &current_particle,
                      size_t iteration_counter)
{
  this->specific_mcmc_adapt(current_particle,
                            iteration_counter);
  current_particle.erase_mcmc_adaptation_info();
}
