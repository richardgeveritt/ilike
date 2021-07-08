#include "smc.h"

//#include "function_pointers.h"
#include "utils.h"
//#include "likelihood_estimator.h"
//#include "likelihood_maker.h"
#include "smc_worker.h"
#include "smc_output.h"

SMC::SMC()
  :LikelihoodEstimator()
{
  //this->output = NULL;
}

SMC::SMC(RandomNumberGenerator* rng_in,
         size_t* seed_in,
         const Data* data_in,
         size_t number_of_particles_in,
         size_t lag_in,
         size_t lag_proposed_in)
  :LikelihoodEstimator(rng_in, seed_in, data_in)
{
  //this->output = new SMCOutput(lag_in,
  //                             lag_proposed_in);
  this->number_of_particles = number_of_particles_in;
  this->lag = lag_in;
  this->lag_proposed = lag_proposed_in;
}

SMC::SMC(const SMC &another)
  :LikelihoodEstimator(another)
{
  this->make_copy(another);
}

SMC::~SMC(void)
{
  if (this->the_worker!=NULL)
    delete this->the_worker;

  //if (this->output!=NULL)
  //  delete this->output;
}

void SMC::operator=(const SMC &another)
{
  if(this == &another)
    return;

  this->make_copy(another);
}

void SMC::make_copy(const SMC &another)
{
  if (this->the_worker!=NULL)
    this->the_worker = another.the_worker->duplicate();

  //if (this->output!=NULL)
  //  this->output = another.output->smc_duplicate();

  this->number_of_particles = another.number_of_particles;
}

Particles SMC::is_step() const
{
  // The way in which this is done is determined by what is set in model_and_algorithm.
  //this->model_and_algorithm;

  // One choice will use an RcppParallel worker.
  return Particles();
}

// Parameters SMC::single_particle_is_step() const
// {
//   //Parameters result = this->model_and_algorithm->simulate_priors->simulate();
//
//   return Parameters();
//
// }


LikelihoodEstimatorOutput* SMC::initial_simulate(const Parameters &parameters)
{
  SMCOutput* output = new SMCOutput(this, this->lag, this->lag_proposed);
  this->smc_update(output);
  return output;
}
