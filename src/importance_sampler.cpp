#include "importance_sampler.h"
#include "smc_worker.h"
//#include "rcppparallel_smc_worker.h"
#include "sequential_smc_worker.h"
#include "smc_output.h"
#include "exact_likelihood_estimator.h"
#include "parameter_particle_simulator.h"

ImportanceSampler::ImportanceSampler()
  :SMC()
{
}

ImportanceSampler::ImportanceSampler(RandomNumberGenerator* rng_in,
                                     size_t* seed_in,
                                     bool parallel_in,
                                     const Data* data_in,
                                     size_t number_of_particles_in,
                                     SimulateDistributionPtr simulate_distribution_in,
                                     EvaluateLogLikelihoodPtr evaluate_log_likelihood_in)
  :SMC(rng_in, seed_in, data_in, number_of_particles_in, 1, 0)
{
   this->model_and_algorithm.likelihood_estimators.resize(0);
   this->model_and_algorithm.likelihood_estimators.reserve(1);
   this->model_and_algorithm.likelihood_estimators.push_back(new ExactLikelihoodEstimator(rng_in,
                                                                      seed_in,
                                                                      data_in,
                                                                      evaluate_log_likelihood_in));
   
   // Need to construct LikelihoodEstimator to read in to this constructor.
   this->model_and_algorithm.particle_simulator = new ParameterParticleSimulator(simulate_distribution_in,
                                                             this->model_and_algorithm.likelihood_estimators);
   
   //if (parallel_in==TRUE)
   //{
   //    this->the_worker = new RcppParallelSMCWorker(this,
   //                                                 this->model_and_algorithm.particle_simulator);
   // }
   // else
   // {
   this->the_worker = new SequentialSMCWorker(this,
                                                 this->model_and_algorithm.particle_simulator);
   // }

}

//Copy constructor for the ImportanceSampler class.
ImportanceSampler::ImportanceSampler(const ImportanceSampler &another)
  :SMC(another)
{
  this->make_copy(another);
}

//Destructor for the ImportanceSampler class.
ImportanceSampler::~ImportanceSampler(void)
{
}

void ImportanceSampler::operator=(const ImportanceSampler &another)
{
  if(this == &another){ //if a==a
    return;
  }

  SMC::operator=(another);
  this->make_copy(another);
}

SMC* ImportanceSampler::smc_duplicate(void)const
{
  return( new ImportanceSampler(*this));
}

LikelihoodEstimator* ImportanceSampler::duplicate(void)const
{
  return( new ImportanceSampler(*this));
}

void ImportanceSampler::make_copy(const ImportanceSampler &another)
{

}

// void ImportanceSampler::smc_step()
// {
// }
//
// void ImportanceSampler::weight_update()
// {
// }

LikelihoodEstimatorOutput* ImportanceSampler::run()
{
   return this->run(Parameters());
}

LikelihoodEstimatorOutput* ImportanceSampler::run(const Parameters &parameters)
{
   return this->initial_simulate(parameters);
}

void ImportanceSampler::smc_update(SMCOutput* current_state)
{
   this->the_worker->simulate();
   Particles particles(Particles(this->the_worker->get_particles()));
   current_state->add_particles(particles);
   //this->output->add_proposed_particles(particles);

   //the_worker->simulate_and_weight();
  /*
   unsigned int number_of_points = algorithm["number_of_points"];

   List observed_data = model["observed_data"];

   // Do the initial importance sampling step.

   // Do the simulation.
   SEXP simulate_proposal_SEXP = algorithm["simulate_proposal"];
   SimulateDistributionPtr simulate_proposal = load_simulate_distribution(simulate_proposal_SEXP);

   std::vector<List> proposed_points;
   proposed_points.reserve(number_of_points);
   for (unsigned int i=0; i<number_of_points; ++i)
   {
   proposed_points.push_back(simulate_proposal());
   }

   LikelihoodEstimator* likelihood_estimator = make_likelihood_estimator(model, algorithm);

   std::vector<List> proposed_auxiliary_variables;
   proposed_auxiliary_variables.reserve(number_of_points);

   for (std::vector<List>::const_iterator i=proposed_points.begin(); i!=proposed_points.end(); ++i)
   {
   proposed_auxiliary_variables.push_back(likelihood_estimator->simulate_auxiliary_variables(*i));
   }

   likelihood_estimator->is_setup_likelihood_estimator(proposed_points,
   proposed_auxiliary_variables);

   arma::colvec log_weights(number_of_points);
   bool prior_is_proposal = algorithm["prior_is_proposal"];
   if (prior_is_proposal==TRUE)
   {
   for (unsigned int i=0; i<number_of_points; ++i)
   {
   log_weights[i] = likelihood_estimator->estimate_log_likelihood(proposed_points[i], proposed_auxiliary_variables[i]);
   }
   }
   else
   {
   SEXP evaluate_log_prior_SEXP = model["evaluate_log_prior"];
   EvaluateLogDistributionPtr evaluate_log_prior = load_evaluate_log_distribution(evaluate_log_prior_SEXP);

   SEXP evaluate_log_proposal_SEXP = algorithm["evaluate_log_proposal"];
   EvaluateLogDistributionPtr evaluate_log_proposal = load_evaluate_log_distribution(evaluate_log_proposal_SEXP);

   for (unsigned int i=0; i<number_of_points; ++i)
   {
   log_weights[i] = likelihood_estimator->estimate_log_likelihood(proposed_points[i], proposed_auxiliary_variables[i]) + evaluate_log_prior(proposed_points[i]) - evaluate_log_proposal(proposed_points[i]);
   }
   }

   if (likelihood_estimator != NULL)
   delete likelihood_estimator;

   return List::create(Named("proposed_points") = proposed_points,
   Named("proposed_auxiliary_variables") = wrap(proposed_auxiliary_variables),
   Named("log_weights") = log_weights,
   Named("log_normalising_constant") = log_sum_exp(log_weights));
   */

}
