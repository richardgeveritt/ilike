#include "smc.h"
#include "progress_reporter.h"
#include <mutex>
#include "chunk_parallel.h"

#include "utils.h"
#include "smc_worker.h"
#include "sequential_smc_worker.h"
#include "smc_output.h"
#include "smc_criterion.h"
#include "ess_smc_criterion.h"
#include "move_output.h"
#include "particle_simulator.h"
#include "smc_termination.h"
#include "ess_smc_criterion.h"
#include "factor_variables.h"
#include "transform.h"
#include "factors.h"

namespace ilike
{
SMC::SMC()
:LikelihoodEstimator()
{
  this->particle_simulator.reset();
  this->setup_default_ancestor_variables();
  this->initialised = false;
  this->the_worker.reset();
  this->resampling_criterion.reset();
}

SMC::SMC(RandomNumberGenerator* rng_in,
         size_t* seed_in,
         Data* data_in,
         const Parameters &algorithm_parameters_in,
         size_t number_of_particles_in,
         size_t lag_in,
         size_t lag_proposed_in,
         const std::vector<const ProposalKernel*> proposals_in,
         double resampling_desired_ess_in,
         bool proposal_is_evaluated_in,
         //EvaluateLogDistributionPtr evaluate_log_proposal_in,
         bool smcfixed_flag_in,
         bool sequencer_limit_is_fixed_in,
         bool transform_proposed_particles,
         const std::string &results_name_in)
:LikelihoodEstimator(rng_in, seed_in, data_in, algorithm_parameters_in, smcfixed_flag_in)
{
  this->number_of_particles = number_of_particles_in;
  this->lag = lag_in;
  this->lag_proposed = lag_proposed_in;
  this->sequencer_limit_is_fixed = sequencer_limit_is_fixed_in;
  this->proposal_is_evaluated = proposal_is_evaluated_in;
  this->resampling_criterion.reset(new ESSSMCCriterion(resampling_desired_ess_in));
  this->proposed_particles_inputted = false;
  this->particle_simulator.reset();
  this->setup_default_ancestor_variables();
  this->results_name = results_name_in;
  
  this->initialised = false;
  // Set up worker?
  
  this->proposals_to_find_gradient_for = proposals_in;
  if (transform_proposed_particles)
    this->proposals_to_transform_for = proposals_in;
  else
    this->proposals_to_transform_for = std::vector<const ProposalKernel*>();
  
  this->the_worker.reset();
}

SMC::SMC(RandomNumberGenerator* rng_in,
         size_t* seed_in,
         Data* data_in,
         const Parameters &algorithm_parameters_in,
         size_t number_of_particles_in,
         size_t lag_in,
         size_t lag_proposed_in,
         const std::vector<const ProposalKernel*> proposals_in,
         SMCCriterion* adaptive_resampling_method_in,
         bool proposal_is_evaluated_in,
         //EvaluateLogDistributionPtr evaluate_log_proposal_in,
         bool smcfixed_flag_in,
         bool sequencer_limit_is_fixed_in,
         bool transform_proposed_particles,
         const std::string &results_name_in)
:LikelihoodEstimator(rng_in, seed_in, data_in, algorithm_parameters_in, smcfixed_flag_in)
{
  this->number_of_particles = number_of_particles_in;
  this->lag = lag_in;
  this->lag_proposed = lag_proposed_in;
  this->sequencer_limit_is_fixed = sequencer_limit_is_fixed_in;
  this->proposal_is_evaluated = proposal_is_evaluated_in;
  this->resampling_criterion.reset(adaptive_resampling_method_in);
  this->proposed_particles_inputted = false;
  this->particle_simulator.reset();
  this->setup_default_ancestor_variables();
  this->results_name = results_name_in;
  
  this->initialised = false;
  // Set up worker?
  
  this->proposals_to_find_gradient_for = proposals_in;
  if (transform_proposed_particles)
    this->proposals_to_transform_for = proposals_in;
  else
    this->proposals_to_transform_for = std::vector<const ProposalKernel*>();
  
  this->the_worker.reset();
}

SMC::SMC(const SMC &another)
:LikelihoodEstimator(another)
{
  this->make_copy(another);
}

// (defined here, not in smc.h, where MoveOutput is only declared: deleting
// through an incomplete type would not run its virtual destructor)
MoveOutput* SMC::move_reusing(RandomNumberGenerator &rng,
                              const Particle &particle,
                              MoveOutput* spare) const
{
  delete spare;
  return this->move(rng, particle);
}

SMC::~SMC()
{
  
  this->the_worker.reset();
  
  this->resampling_criterion.reset();
  
  this->particle_simulator.reset();

}

void SMC::operator=(const SMC &another)
{
  if(this == &another)
    return;
  
  this->the_worker.reset();
  
  this->resampling_criterion.reset();
  
  this->particle_simulator.reset();

  this->make_copy(another);
}

void SMC::make_copy(const SMC &another)
{
  this->write_factor_output = another.write_factor_output;
  this->compact_stored_output = another.compact_stored_output;
  this->chunk_size = another.chunk_size;
  this->force_rich_generations = another.force_rich_generations;
  this->parallel = another.parallel;
  this->number_of_threads = another.number_of_threads;
  this->waste_free_chains = another.waste_free_chains;

  if (another.the_worker!=NULL)
    this->the_worker.reset(another.the_worker->duplicate());
  else
    this->the_worker.reset();
  
  if (another.resampling_criterion!=NULL)
    this->resampling_criterion.reset(another.resampling_criterion->duplicate());
  else
    this->resampling_criterion.reset();
  
  if (another.particle_simulator!=NULL)
    this->particle_simulator.reset(another.particle_simulator->duplicate());
  else
    this->particle_simulator.reset();

  this->sequencer = another.sequencer;
  // B29: the copy's worker and sequencer work for the copy, and its schedule
  // variables are its own (Sequencer's copy shares the matrices, and the
  // sequencer writes into them in place)
  if (this->the_worker!=NULL)
    this->the_worker->set_smc(this);
  this->sequencer.set_worker(this->the_worker.get());
  this->sequencer.schedule_parameters = another.sequencer.schedule_parameters.deep_copy();
  
  this->sequencer_limit_is_fixed = another.sequencer_limit_is_fixed;
  this->number_of_particles = another.number_of_particles;
  this->lag = another.lag;
  this->lag_proposed = another.lag_proposed;
  this->proposal_is_evaluated = another.proposal_is_evaluated;
  this->proposed_particles_inputted = another.proposed_particles_inputted;
  this->initial_particles = another.initial_particles;
  this->vector_variables = another.vector_variables;
  this->any_variables = another.any_variables;
  this->vector_variable_sizes = another.vector_variable_sizes;
  this->default_ancestor_variables = another.default_ancestor_variables;
  this->log_probabilities_of_initial_values = another.log_probabilities_of_initial_values;
  this->initialised = another.initialised;
  this->results_name = another.results_name;
  
  this->proposals_to_transform_for = another.proposals_to_transform_for;
  this->proposals_to_find_gradient_for = another.proposals_to_find_gradient_for;
}

//Particles SMC::is_step() const
//{
// The way in which this is done is determined by what is set in model_and_algorithm.
//this->model_and_algorithm;

// One choice will use an RcppParallel worker.
//  return Particles();
//}

// Parameters SMC::single_particle_is_step() const
// {
//   //Parameters result = this->model_and_algorithm->simulate_priors->simulate();
//
//   return Parameters();
//
// }

SMCOutput* SMC::run()
{
  
  if (this->initialised==false)
  {
    this->setup();
    this->initialised = true;
  }
  
  return this->specific_run();
}

SMCOutput* SMC::run(const Parameters &conditioned_on_parameters)
{
  this->setup(conditioned_on_parameters);
  return this->specific_run(conditioned_on_parameters);
}

LikelihoodEstimatorOutput* SMC::initialise()
{
  return this->initialise_smc();
}

SMCOutput* SMC::initialise_smc()
{
  if (this->initialised==false)
  {
    this->setup();
    this->initialised = true;
  }
  
  return this->specific_initialise_smc();
}

void SMC::setup()
{
  this->setup_variables();
}

void SMC::setup(const Parameters &parameters)
{
  this->setup_variables(parameters);
}

void SMC::setup_variables()
{
  Parameters dummy_parameters;
  
  if (this->proposed_particles_inputted)
  {
    dummy_parameters = this->initial_particles[0];
    dummy_parameters.merge_with_fixed(this->sequencer.schedule_parameters);
  }
  else
  {
    dummy_parameters = std::move(this->particle_simulator->simulate(*this->rng,
                                                                    this->factors.get(),
                                                                    &this->proposals_to_transform_for,
                                                                    &this->proposals_to_find_gradient_for,
                                                                    this->sequencer.schedule_parameters).parameters);
  }
  
  this->vector_variables = dummy_parameters.get_nonfixed_vector_variables();
  this->any_variables = dummy_parameters.get_nonfixed_any_variables();
  
  this->vector_variable_sizes = dummy_parameters.get_variable_n_elems(this->vector_variables);
  
  if (this->factors!=NULL)
    this->factors->setup(dummy_parameters);
  
}

void SMC::setup_variables(const Parameters &parameters)
{
  Parameters dummy_parameters;
  if (this->proposed_particles_inputted)
  {
    dummy_parameters = this->initial_particles[0];
    dummy_parameters.merge_with_fixed(parameters);
    dummy_parameters.merge_with_fixed(this->sequencer.schedule_parameters);
  }
  else
  {
    dummy_parameters = std::move(this->particle_simulator->simulate(*this->rng,
                                                                    this->factors.get(),
                                                                    &this->proposals_to_transform_for,
                                                                    &this->proposals_to_find_gradient_for,
                                                                    parameters,
                                                                    this->sequencer.schedule_parameters).parameters);
  }
  this->vector_variables = dummy_parameters.get_nonfixed_vector_variables();
  this->any_variables = dummy_parameters.get_nonfixed_any_variables();
  
  this->vector_variable_sizes = dummy_parameters.get_variable_n_elems(this->vector_variables);
  
  if (this->factors!=NULL)
    this->factors->setup(dummy_parameters);
}

void SMC::simulate_generation_in_chunks(SMCOutput *current_state,
                                       const Parameters* conditioned_on_parameters,
                                       bool summarise,
                                       const Index* index,
                                       bool keep_points)
  {
    // What simulate_smc + evaluate_smcfixed_part_smc + the weighting step of
    // the bisection do (SequentialSMCWorker::specific_simulate,
    // SMCWorker::simulate, smcfixed_weight, smcadaptive_given_smcfixed_weight),
    // B points at a time. Per chunk the phases keep their order (simulate all
    // B, then the smcfixed parts, then the weights or summaries). Chunk k uses
    // its own random-number stream (chunk_streams). Three modes (dev/plan.md
    // 2.6a G4, 2.6b G5):
    //  - independent: no factor uses random numbers or keeps run state; every
    //    chunk (chunk 0 included) may run in parallel, sharing the factors;
    //  - copies: some factors run an inner algorithm; each group of chunks
    //    gets its own copy of the factors, drawing from its own generator, and
    //    each point's inner runs draw from a seed counter set from the point's
    //    chunk stream (and kept with the point between its simulation and its
    //    smcfixed part), so the results do not depend on the threads;
    //  - shared: anything else; sequential, in the original order (the seed
    //    bookkeeping of SMCWorker::simulate after chunk 0's simulations), so
    //    that with one chunk the results are those of the per-particle path.
    // In the first two modes the seed bookkeeping is done first (after the
    // chunk streams have been taken from the seed: it does not touch them),
    // and the first point to be kept fixes the layout of what is kept.
    SMCWorker* worker = this->the_worker.get();
    size_t N = this->number_of_particles;
    size_t B = ((this->chunk_size == 0) || (this->chunk_size > N)) ? N : this->chunk_size;
    size_t K = (N + B - 1)/B;

    enum Mode { independent, copies, shared };
    Mode mode = shared;
    if ((this->factors != NULL) && this->factors->supports_parallel_chunks())
      mode = independent;
    else if ((this->factors != NULL) && this->factors->supports_parallel_copies())
      mode = copies;
    bool in_parallel = this->parallel && (mode != shared);
    if (this->parallel && !in_parallel)
      this->parallel_refused = true;

    Particles* generation = current_state->add_particles();
    generation->flat = true;
    generation->flat_size = N;
    if (!summarise)
      generation->flat_incremental_log_weights.set_size(N);

    // What is kept of each point: its variables and factor outputs only if
    // keep_points (an estimator inside a factor stores no generations; an SMC
    // keeps its first generation's points to move them).
    bool keep_factor_output = keep_points && this->write_factor_output;
    bool compact_records = false;
    size_t number_of_records = 0;
    int summary_size = -1;
    std::once_flag store_setup;
    std::once_flag summary_setup;

    RandomNumberGenerator first_stream(*worker->get_rng());
    first_stream.seed(worker->get_seed());
    std::vector<RandomNumberGenerator> streams = chunk_streams(first_stream, K);

    auto seed_bookkeeping = [&]()
    {
      // the seed bookkeeping and resampling variables of SMCWorker::simulate
      worker->set_seed(worker->get_seed() + N);
      this->rng->seed(worker->get_seed());
      generation->simulate_resampling_variables(*this->rng);
      worker->set_seed(worker->get_seed() + N);
    };

    // What a group of chunks works with: the factors (the main ones, or the
    // group's copy), its particle objects, and in copies mode the generator
    // and seed counter its copy draws from.
    struct Workspace
    {
      std::unique_ptr<Factors> own_factors;
      const Factors* factors = NULL;
      RandomNumberGenerator rng;
      size_t seed = 0;
      std::vector<size_t> point_seeds;
      std::vector<Particle> slots;
      std::vector<double> records;
    };
    auto make_workspace = [&](Workspace &w)
    {
      if (mode == copies)
      {
        w.own_factors.reset(this->factors->duplicate());
        w.own_factors->redirect_random_numbers(&w.rng, &w.seed);
        w.factors = w.own_factors.get();
        w.point_seeds.assign(B, 0);
      }
      else
      {
        w.factors = this->factors.get();
      }
      w.slots.resize(B);
    };

    // One chunk: simulate, evaluate and keep points [first, first + n).
    auto run_chunk = [&](size_t k, Workspace &w)
    {
      size_t first = k*B;
      size_t n = std::min(B, N - first);
      RandomNumberGenerator local_rng(streams[k]);
      std::vector<Particle> &slots = w.slots;

      for (size_t j = 0; j < n; ++j)
      {
        if (mode == copies)
          w.seed = size_t(local_rng());
        if (conditioned_on_parameters == NULL)
          this->particle_simulator->simulate_into(local_rng,
                                                  slots[j],
                                                  w.factors,
                                                  &this->proposals_to_transform_for,
                                                  &this->proposals_to_find_gradient_for,
                                                  this->sequencer.schedule_parameters);
        else
          this->particle_simulator->simulate_into(local_rng,
                                                  slots[j],
                                                  w.factors,
                                                  &this->proposals_to_transform_for,
                                                  &this->proposals_to_find_gradient_for,
                                                  *conditioned_on_parameters,
                                                  this->sequencer.schedule_parameters);
        if (this->proposal_is_evaluated == true)
          slots[j].previous_target_evaluated = this->particle_simulator->evaluate(slots[j]);
        if (mode == copies)
          w.point_seeds[j] = w.seed;
      }

      if ((mode == shared) && (k == 0))
        seed_bookkeeping();

      for (size_t j = 0; j < n; ++j)
      {
        if (mode == copies)
          w.seed = w.point_seeds[j];
        slots[j].evaluate_smcfixed_part_of_likelihoods(index);
      }

      if (summarise)
      {
        for (size_t j = 0; j < n; ++j)
        {
          const FactorVariables* factor_variables = slots[j].factor_variables.get();
          int size = (factor_variables == NULL) ? -1 : factor_variables->adaptive_summary_size(index);
          std::call_once(summary_setup, [&]()
          {
            if (size < 0)
              Rcpp::stop("SMC::simulate_generation_in_chunks: a factor has no adaptive summary.");
            summary_size = size;
            generation->adaptive_summaries.set_size(size_t(summary_size), N);
            generation->adaptive_previous_targets.set_size(N);
            generation->flat_summary_evaluator.reset(factor_variables->duplicate());
            generation->flat_summary_evaluator->set_particle(NULL);
            if ((mode == copies) && !generation->flat_summary_evaluator->repoint_to(this->factors.get()))
              Rcpp::stop("SMC::simulate_generation_in_chunks: the factors could not be repointed.");
            generation->adaptive_summary_evaluator = generation->flat_summary_evaluator.get();
            generation->adaptive_summary_index = index;
            generation->adaptive_summary_state = 1;
          });
          if (size != summary_size)
            Rcpp::stop("SMC::simulate_generation_in_chunks: the points' adaptive summaries differ in size.");
          if (summary_size > 0)
            factor_variables->adaptive_summary(index, generation->adaptive_summaries.colptr(first + j));
          generation->adaptive_previous_targets[first + j] = slots[j].previous_target_evaluated;
        }
      }
      else
      {
        for (size_t j = 0; j < n; ++j)
        {
          double a = slots[j].evaluate_smcadaptive_part_given_smcfixed_likelihoods(index);
          double b = slots[j].previous_target_evaluated;
          generation->flat_incremental_log_weights[first + j] = (a == -arma::datum::inf) ? -arma::datum::inf : a - b;
        }
      }

      if (!keep_points)
        return;
      for (size_t j = 0; j < n; ++j)
      {
        Particle &slot = slots[j];
        size_t point = first + j;
        std::call_once(store_setup, [&]()
        {
          if (keep_factor_output && (slot.factor_variables != NULL))
          {
            compact_records = slot.factor_variables->has_compact_values();
            if (compact_records)
            {
              number_of_records = slot.factor_variables->number_of_compact_values();
              generation->flat_points.set_number_of_extra_columns(number_of_records);
              generation->flat_record_writer.reset(slot.factor_variables->duplicate());
              generation->flat_record_writer->set_particle(NULL);
              // summarised points have not been weighted yet: their records
              // are recomputed once the target has been chosen
              generation->flat_records_from_summaries = summarise;
              if ((mode == copies) && !generation->flat_record_writer->repoint_to(this->factors.get()))
                Rcpp::stop("SMC::simulate_generation_in_chunks: the factors could not be repointed.");
            }
            else
            {
              generation->flat_factor_variables.clear();
              generation->flat_factor_variables.resize(N);
            }
          }
          generation->flat_points.prepare(slot.parameters, N);
        });
        if (compact_records)
        {
          w.records.resize(number_of_records);
          slot.factor_variables->compact_values(w.records.data());
          generation->flat_points.set_point(point, slot.parameters, w.records.data());
        }
        else
        {
          generation->flat_points.set_point(point, slot.parameters);
          if (keep_factor_output && (slot.factor_variables != NULL))
          {
            // keep this point's factor outputs whole (repointed to the main
            // factors, which outlive the group's copy); the slot gets new ones
            slot.factor_variables->set_particle(NULL);
            if ((mode == copies) && !slot.factor_variables->repoint_to(this->factors.get()))
              Rcpp::stop("SMC::simulate_generation_in_chunks: the factors could not be repointed.");
            generation->flat_factor_variables[point] = std::move(slot.factor_variables);
          }
        }
      }
    };

    if (mode == shared)
    {
      Workspace w;
      make_workspace(w);
      for (size_t k = 0; k < K; ++k)
      {
        run_chunk(k, w);
        if (this->progress_counts_points)
          ProgressReporter::instance().update(double(std::min(N, (k + 1)*B)));
        check_interrupt();
      }
    }
    else
    {
      seed_bookkeeping();
      std::function<void(size_t)> after_wave;
      if (this->progress_counts_points)
        after_wave = [&](size_t chunks_done) { ProgressReporter::instance().update(double(std::min(N, chunks_done*B))); };
      run_chunk_groups(0, K, in_parallel, this->number_of_threads,
                       [&](size_t group, size_t first_chunk, size_t end_chunk)
                       {
                         Workspace w;
                         make_workspace(w);
                         for (size_t k = first_chunk; k < end_chunk; ++k)
                           run_chunk(k, w);
                         // the slots' factor variables belong to the group's
                         // copy of the factors: delete them first
                         w.slots.clear();
                       },
                       after_wave);
    }

    generation->initialise_weights();
  }

void SMC::weight_from_summaries_if_possible(const Index* index,
                                            Particles &current_particles)
{
  if (this->force_rich_generations ||
      (dynamic_cast<SequentialSMCWorker*>(this->the_worker.get()) == NULL) ||
      !current_particles.prepare_adaptive_summaries(index))
  {
    this->the_worker->smcadaptive_given_smcfixed_weight(index, current_particles);
    return;
  }
  // as SequentialSMCWorker::smcadaptive_given_smcfixed_weight
  size_t n = current_particles.size();
  arma::colvec weights(n);
  const FactorVariables* evaluator = current_particles.adaptive_summary_evaluator;
  const Parameters &schedule_parameters = this->sequencer.schedule_parameters;
  std::vector<double> targets(n);
  // in parallel when asked: each particle's value depends only on its own
  // summary (no random numbers, nothing shared is written), so the values
  // are the same in any order
  size_t B = ((this->chunk_size == 0) || (this->chunk_size > n)) ? n : this->chunk_size;
  size_t K = (B == 0) ? 0 : (n + B - 1)/B;
  size_t stride = current_particles.adaptive_summaries.n_rows;
  const double* summaries = current_particles.adaptive_summaries.memptr();
  run_chunk_groups(0, K, this->parallel, this->number_of_threads,
                   [&](size_t group, size_t first_chunk, size_t end_chunk)
                   {
                     size_t first = first_chunk*B;
                     size_t last = std::min(n, end_chunk*B);
                     if (last > first)
                       evaluator->smcadaptive_from_summaries(index,
                                                             summaries + first*stride,
                                                             stride,
                                                             last - first,
                                                             schedule_parameters,
                                                             targets.data() + first);
                   });
  for (size_t i = 0; i < n; ++i)
  {
    double a = targets[i];
    double b = current_particles.adaptive_previous_targets[i];
    if (a==-arma::datum::inf)
      weights[i] = -arma::datum::inf;
    else
      weights[i] = a - b;
  }
  this->the_worker->set_unnormalised_log_incremental_weights(weights);
  current_particles.adaptive_summaries_used = true;
  current_particles.last_trial_schedule = schedule_parameters.deep_copy();
}

// Write the values of `values` into the matrices of `target` (in place: the
// particles share the sequencer's schedule-variable matrices).
static void assign_schedule_values(Parameters &target, const Parameters &values)
{
  std::vector<std::string> names = values.get_vector_variables();
  for (auto i = names.begin(); i != names.end(); ++i)
    target[*i] = values[*i];
}

void SMC::finish_weighting_from_summaries(const Index* index,
                                          Particles &current_particles)
{
  if (!current_particles.adaptive_summaries_used)
    return;
  // a flat generation has no cached values to bring up to date; it keeps
  // last_trial_schedule, at which its next moves restore its factor values;
  // records taken before its target was chosen are recomputed as the values
  // the last weighting would have left (factors restored from the summaries)
  if (current_particles.is_flat())
  {
    if (current_particles.flat_records_from_summaries &&
        (current_particles.flat_record_writer != NULL) &&
        (current_particles.adaptive_summary_state == 1))
    {
      FactorVariables* writer = current_particles.flat_record_writer.get();
      std::vector<double> records(writer->number_of_compact_values());
      for (size_t i = 0; i < current_particles.size(); ++i)
      {
        if (!writer->restore_from_summary(index,
                                          NULL,
                                          current_particles.adaptive_summaries.colptr(i),
                                          current_particles.last_trial_schedule))
          Rcpp::stop("SMC::finish_weighting_from_summaries: the factor records could not be recomputed.");
        writer->compact_values(records.data());
        current_particles.flat_points.set_extras(i, records.data());
      }
    }
    return;
  }
  current_particles.adaptive_summaries_used = false;
  Parameters chosen = this->sequencer.schedule_parameters.deep_copy();
  assign_schedule_values(this->sequencer.schedule_parameters, current_particles.last_trial_schedule);
  this->the_worker->smcadaptive_given_smcfixed_weight(index, current_particles);
  current_particles.update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
  assign_schedule_values(this->sequencer.schedule_parameters, chosen);
}

void SMC::simulate_proposal(SMCOutput* current_state)
{
  // Simulate from proposal.
  Particles* next_particles = current_state->add_particles();
  
  if (!this->proposed_particles_inputted)
  {
    this->the_worker->simulate(next_particles);
  }
  else
  {
    next_particles->setup(this->initial_particles,
                          this->log_probabilities_of_initial_values,
                          this->factors.get(),
                          &this->proposals_to_transform_for,
                          &this->proposals_to_find_gradient_for,
                          this->sequencer.schedule_parameters);
  }
  
  current_state->back().initialise_weights();
}

LikelihoodEstimatorOutput* SMC::initialise(const Parameters &conditioned_on_parameters)
{
  return this->initialise_smc(conditioned_on_parameters);
}

SMCOutput* SMC::initialise_smc(const Parameters &conditioned_on_parameters)
{
  
  if (this->initialised==false)
  {
    this->setup(conditioned_on_parameters);
    this->initialised = true;
  }
  
  return this->specific_initialise_smc(conditioned_on_parameters);
}

void SMC::simulate_proposal(SMCOutput* current_state,
                            const Parameters &conditioned_on_parameters)
{
  // Simulate from proposal.
  Particles* next_particles = current_state->add_particles();
  
  if (!this->proposed_particles_inputted)
  {
    this->the_worker->simulate(next_particles,
                               conditioned_on_parameters);
  }
  else
  {
    next_particles->setup(this->initial_particles,
                          this->log_probabilities_of_initial_values,
                          this->factors.get(),
                          &this->proposals_to_transform_for,
                          &this->proposals_to_find_gradient_for,
                          conditioned_on_parameters,
                          this->sequencer.schedule_parameters);
  }
  
  current_state->back().initialise_weights();
}

void SMC::resample(SMCOutput* current_state)
{
  auto current_particles_iterator = current_state->end()-1;
  
  // Check criterion to see if we resample.
  if ((*this->resampling_criterion.get())(*current_particles_iterator)<0.0)
  {
    // Sample ancestor variables.
    current_particles_iterator->resampled_flag = true;
    current_particles_iterator->resample();
  }
  else
  {
    current_particles_iterator->ancestor_variables = this->default_ancestor_variables;
  }
}

void SMC::setup_default_ancestor_variables()
{
  // Set ancestor variables to be 1:n.
  this->default_ancestor_variables.clear();
  this->default_ancestor_variables.reserve(this->number_of_particles);
  for (size_t i=0; i<this->number_of_particles; ++i)
  {
    this->default_ancestor_variables.push_back(i);
  }
}

void SMC::specific_change_data(Data* new_data)
{
  
}
}
