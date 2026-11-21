#include "importance_sampler.h"
#include "smc_worker.h"
#include "sequential_smc_worker.h"
#include "smc_output.h"
#include "utils.h"
#include "exact_likelihood_estimator.h"
#include "parameter_particle_simulator.h"
#include "move_output.h"
#include "single_point_move_output.h"
#include "independent_proposal_kernel.h"
#include "custom_distribution_proposal_kernel.h"
#include "vector_factors.h"
#include "vector_index.h"
#include "positive_smc_criterion.h"
#include "factor_variables.h"
#include "chunk_parallel.h"
#include "progress_reporter.h"
#include <mutex>
#include "particle_simulator.h"
#include <memory>

namespace ilike
{
  ImportanceSampler::ImportanceSampler()
      : SMC()
  {
    this->index.reset();
  }

  ImportanceSampler::ImportanceSampler(RandomNumberGenerator *rng_in,
                                       size_t *seed_in,
                                       Data *data_in,
                                       const Parameters &algorithm_parameters_in,
                                       size_t number_of_particles_in,
                                       EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
                                       SimulateDistributionPtr simulate_prior_in,
                                       bool store_output,
                                       bool smcfixed_flag_in,
                                       bool transform_proposed_particles,
                                       bool parallel_in,
                                       size_t grain_size_in,
                                       const std::string &results_name_in)
      : SMC(rng_in, seed_in, data_in, algorithm_parameters_in, number_of_particles_in, size_t(store_output), 0, std::vector<const ProposalKernel *>(), double(number_of_particles_in), false, smcfixed_flag_in, true, transform_proposed_particles, results_name_in)
  {
    this->chunk_size = grain_size_in;
    IndependentProposalKernel *proposal = new CustomDistributionProposalKernel(simulate_prior_in);
    proposal->set_proposal_parameters(&this->algorithm_parameters);

    std::vector<size_t> indices;
    std::vector<LikelihoodEstimator *> likelihood_estimators;
    likelihood_estimators.reserve(1);
    likelihood_estimators.push_back(new ExactLikelihoodEstimator(rng_in,
                                                                 seed_in,
                                                                 data_in,
                                                                 evaluate_log_likelihood_in,
                                                                 true));
    indices.push_back(0);
    this->index.reset(new VectorIndex(indices));

    this->factors.reset(new VectorFactors(likelihood_estimators));

    // Need to construct LikelihoodEstimator to read in to this constructor.
    this->particle_simulator.reset(new ParameterParticleSimulator(proposal,
                                                              likelihood_estimators));

    this->parallel = parallel_in;
    this->the_worker.reset(new SequentialSMCWorker(this));

    std::vector<double> schedule_in;
    schedule_in.push_back(0.0);
    schedule_in.push_back(1.0);
    SMCCriterion *smc_criterion = new PositiveSMCCriterion();
    this->sequencer = Sequencer(this->the_worker.get(),
                                schedule_in,
                                "",
                                25,
                                smc_criterion);
  }

  ImportanceSampler::ImportanceSampler(RandomNumberGenerator *rng_in,
                                       size_t *seed_in,
                                       Data *data_in,
                                       const Parameters &algorithm_parameters_in,
                                       size_t number_of_particles_in,
                                       EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
                                       EvaluateLogDistributionPtr evaluate_log_prior_in,
                                       SimulateDistributionPtr simulate_proposal_in,
                                       EvaluateLogDistributionPtr evaluate_log_proposal_in,
                                       bool store_output,
                                       bool smcfixed_flag_in,
                                       bool transform_proposed_particles,
                                       bool parallel_in,
                                       size_t grain_size_in,
                                       const std::string &results_name_in)
      : SMC(rng_in, seed_in, data_in, algorithm_parameters_in, number_of_particles_in, size_t(store_output), 0, std::vector<const ProposalKernel *>(), double(number_of_particles_in), true, smcfixed_flag_in, true, transform_proposed_particles, results_name_in)
  {
    this->chunk_size = grain_size_in;
    IndependentProposalKernel *proposal = new CustomDistributionProposalKernel(simulate_proposal_in,
                                                                               evaluate_log_proposal_in);
    proposal->set_proposal_parameters(&this->algorithm_parameters);

    std::vector<LikelihoodEstimator *> likelihood_estimators;
    std::vector<size_t> indices;
    likelihood_estimators.reserve(1);
    likelihood_estimators.push_back(new ExactLikelihoodEstimator(rng_in,
                                                                 seed_in,
                                                                 data_in,
                                                                 evaluate_log_prior_in,
                                                                 evaluate_log_likelihood_in,
                                                                 true));
    indices.push_back(0);
    this->index.reset(new VectorIndex(indices));

    this->factors.reset(new VectorFactors(likelihood_estimators));

    // Need to construct LikelihoodEstimator to read in to this constructor.
    this->particle_simulator.reset(new ParameterParticleSimulator(proposal,
                                                              likelihood_estimators));

    this->parallel = parallel_in;
    this->the_worker.reset(new SequentialSMCWorker(this));

    std::vector<double> schedule_in;
    schedule_in.push_back(0.0);
    schedule_in.push_back(1.0);
    SMCCriterion *smc_criterion = new PositiveSMCCriterion();
    this->sequencer = Sequencer(this->the_worker.get(),
                                schedule_in,
                                "",
                                25,
                                smc_criterion);
  }

  ImportanceSampler::ImportanceSampler(RandomNumberGenerator *rng_in,
                                       size_t *seed_in,
                                       Data *data_in,
                                       const Parameters &algorithm_parameters_in,
                                       size_t number_of_particles_in,
                                       EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
                                       EvaluateLogDistributionPtr evaluate_log_prior_in,
                                       IndependentProposalKernel *proposal_in,
                                       bool store_output,
                                       bool smcfixed_flag_in,
                                       bool transform_proposed_particles,
                                       bool parallel_in,
                                       size_t grain_size_in,
                                       const std::string &results_name_in)
      : SMC(rng_in, seed_in, data_in, algorithm_parameters_in, number_of_particles_in, size_t(store_output), 0, std::vector<const ProposalKernel *>(), double(number_of_particles_in), true, smcfixed_flag_in, true, transform_proposed_particles, results_name_in)
  {
    this->chunk_size = grain_size_in;
    proposal_in->set_proposal_parameters(&this->algorithm_parameters);
    std::vector<size_t> indices;
    std::vector<LikelihoodEstimator *> likelihood_estimators;
    likelihood_estimators.reserve(1);
    likelihood_estimators.push_back(new ExactLikelihoodEstimator(rng_in,
                                                                 seed_in,
                                                                 data_in,
                                                                 evaluate_log_prior_in,
                                                                 evaluate_log_likelihood_in,
                                                                 true));

    indices.push_back(0);
    this->index.reset(new VectorIndex(indices));

    this->factors.reset(new VectorFactors(likelihood_estimators));

    // Need to construct LikelihoodEstimator to read in to this constructor.
    this->particle_simulator.reset(new ParameterParticleSimulator(proposal_in,
                                                              likelihood_estimators));

    this->parallel = parallel_in;
    this->the_worker.reset(new SequentialSMCWorker(this));

    std::vector<double> schedule_in;
    schedule_in.push_back(0.0);
    schedule_in.push_back(1.0);
    SMCCriterion *smc_criterion = new PositiveSMCCriterion();
    this->sequencer = Sequencer(this->the_worker.get(),
                                schedule_in,
                                "",
                                25,
                                smc_criterion);
  }

  ImportanceSampler::ImportanceSampler(RandomNumberGenerator *rng_in,
                                       size_t *seed_in,
                                       Data *data_in,
                                       const Parameters &algorithm_parameters_in,
                                       size_t number_of_particles_in,
                                       const std::string &target_variable_in,
                                       const std::vector<LikelihoodEstimator *> &likelihood_estimators_in,
                                       IndependentProposalKernel *proposal_in,
                                       Index* without_cancelled_index,
                                       bool proposal_is_evaluated_in,
                                       bool store_output,
                                       bool smcfixed_flag_in,
                                       bool sequencer_limit_is_fixed_in,
                                       bool transform_proposed_particles,
                                       bool parallel_in,
                                       size_t grain_size_in,
                                       const std::string &results_name_in)
      : SMC(rng_in, seed_in, data_in, algorithm_parameters_in, number_of_particles_in, size_t(store_output), 0, std::vector<const ProposalKernel *>(), double(number_of_particles_in), proposal_is_evaluated_in, smcfixed_flag_in, sequencer_limit_is_fixed_in, transform_proposed_particles, results_name_in)
  {
    this->chunk_size = grain_size_in;
    proposal_in->set_proposal_parameters(&this->algorithm_parameters);

    this->index.reset(without_cancelled_index);

    this->factors.reset(new VectorFactors(likelihood_estimators_in));

    // Need to construct LikelihoodEstimator to read in to this constructor.
    this->particle_simulator.reset(new ParameterParticleSimulator(proposal_in,
                                                              likelihood_estimators_in));

    this->parallel = parallel_in;
    this->the_worker.reset(new SequentialSMCWorker(this));

    std::vector<double> schedule_in;
    schedule_in.push_back(1.0);
    SMCCriterion *smc_criterion = new PositiveSMCCriterion();
    this->sequencer = Sequencer(this->the_worker.get(),
                                schedule_in,
                                target_variable_in,
                                100,
                                smc_criterion);
  }

  ImportanceSampler::ImportanceSampler(RandomNumberGenerator *rng_in,
                                       size_t *seed_in,
                                       Data *data_in,
                                       const Parameters &algorithm_parameters_in,
                                       size_t number_of_particles_in,
                                       const std::string &target_variable_in,
                                       double target_value_in,
                                       const std::vector<LikelihoodEstimator *> &likelihood_estimators_in,
                                       IndependentProposalKernel *proposal_in,
                                       Index* without_cancelled_index,
                                       bool proposal_is_evaluated_in,
                                       bool store_output,
                                       bool smcfixed_flag_in,
                                       bool sequencer_limit_is_fixed_in,
                                       bool transform_proposed_particles,
                                       bool parallel_in,
                                       size_t grain_size_in,
                                       const std::string &results_name_in)
      : SMC(rng_in, seed_in, data_in, algorithm_parameters_in, number_of_particles_in, size_t(store_output), 0, std::vector<const ProposalKernel *>(), double(number_of_particles_in), proposal_is_evaluated_in, smcfixed_flag_in, sequencer_limit_is_fixed_in, transform_proposed_particles, results_name_in)
  {
    this->chunk_size = grain_size_in;
    proposal_in->set_proposal_parameters(&this->algorithm_parameters);

    this->index.reset(without_cancelled_index);

    this->factors.reset(new VectorFactors(likelihood_estimators_in));

    // Need to construct LikelihoodEstimator to read in to this constructor.
    this->particle_simulator.reset(new ParameterParticleSimulator(proposal_in,
                                                              likelihood_estimators_in));

    this->parallel = parallel_in;
    this->the_worker.reset(new SequentialSMCWorker(this));

    std::vector<double> schedule_in;
    schedule_in.push_back(target_value_in);
    SMCCriterion *smc_criterion = new PositiveSMCCriterion();
    this->sequencer = Sequencer(this->the_worker.get(),
                                schedule_in,
                                target_variable_in,
                                25,
                                smc_criterion);
  }

  ImportanceSampler::ImportanceSampler(RandomNumberGenerator *rng_in,
                                       size_t *seed_in,
                                       Data *data_in,
                                       const Parameters &algorithm_parameters_in,
                                       size_t number_of_particles_in,
                                       const std::string &target_variable_in,
                                       const std::vector<std::string> &measurement_variables_in,
                                       const std::string &scale_variable_in,
                                       const arma::colvec &scale_in,
                                       const std::vector<LikelihoodEstimator *> &likelihood_estimators_in,
                                       IndependentProposalKernel *proposal_in,
                                       bool proposal_is_evaluated_in,
                                       bool store_output,
                                       bool smcfixed_flag_in,
                                       bool sequencer_limit_is_fixed_in,
                                       bool transform_proposed_particles,
                                       bool parallel_in,
                                       size_t grain_size_in,
                                       const std::string &results_name_in)
      : SMC(rng_in, seed_in, data_in, algorithm_parameters_in, number_of_particles_in, size_t(store_output), 0, std::vector<const ProposalKernel *>(), double(number_of_particles_in), proposal_is_evaluated_in, smcfixed_flag_in, sequencer_limit_is_fixed_in, transform_proposed_particles, results_name_in)
  {
    this->chunk_size = grain_size_in;
    proposal_in->set_proposal_parameters(&this->algorithm_parameters);

    std::vector<size_t> indices;
    indices.reserve(likelihood_estimators_in.size());
    for (size_t i = 0; i < likelihood_estimators_in.size(); ++i)
      indices.push_back(i);
    this->index.reset(new VectorIndex(indices));

    this->factors.reset(new VectorFactors(likelihood_estimators_in));

    // Need to construct LikelihoodEstimator to read in to this constructor.
    this->particle_simulator.reset(new ParameterParticleSimulator(proposal_in,
                                                              likelihood_estimators_in));

    this->parallel = parallel_in;
    this->the_worker.reset(new SequentialSMCWorker(this));

    std::vector<double> schedule_in;
    schedule_in.push_back(1.0);
    SMCCriterion *smc_criterion = new PositiveSMCCriterion();
    this->sequencer = Sequencer(this->the_worker.get(),
                                schedule_in,
                                target_variable_in,
                                100,
                                smc_criterion);

    // sort out ABC scaling

    size_t data_size = 0;
    for (auto i = measurement_variables_in.begin();
         i != measurement_variables_in.end();
         ++i)
    {
      data_size = data_size + (*this->data)[*i].n_elem;
    }

    // if scale variable is not set, we are not using a scale at all, and we set it to be ones
    if (scale_variable_in == "")
    {
      this->sequencer.scale_variable = "default_scale";
      arma::colvec scale = arma::colvec(data_size);
      scale.fill(1.0);
      this->sequencer.schedule_parameters[this->sequencer.scale_variable] = scale;
      this->sequencer.find_scale = false;
    }
    else // if the scale variable is set, then either we have read in a valid scale, or it must be estimated through simulation
    {
      this->sequencer.scale_variable = scale_variable_in;

      if (scale_in.n_elem == data_size)
      {
        this->sequencer.schedule_parameters[this->sequencer.scale_variable] = scale_in;
        this->sequencer.find_scale = false;
      }
      else
      {
        Rcpp::stop("Scale variable is set, but the scale is not the correct size (it does not match the dimensions of the observed data).");
      }
    }
  }

  ImportanceSampler::ImportanceSampler(RandomNumberGenerator *rng_in,
                                       size_t *seed_in,
                                       Data *data_in,
                                       const Parameters &algorithm_parameters_in,
                                       size_t number_of_particles_in,
                                       const std::string &target_variable_in,
                                       double target_value_in,
                                       const std::vector<std::string> &measurement_variables_in,
                                       const std::string &scale_variable_in,
                                       const arma::colvec &scale_in,
                                       const std::vector<LikelihoodEstimator *> &likelihood_estimators_in,
                                       IndependentProposalKernel *proposal_in,
                                       bool proposal_is_evaluated_in,
                                       bool store_output,
                                       bool smcfixed_flag_in,
                                       bool sequencer_limit_is_fixed_in,
                                       bool transform_proposed_particles,
                                       bool parallel_in,
                                       size_t grain_size_in,
                                       const std::string &results_name_in)
      : SMC(rng_in, seed_in, data_in, algorithm_parameters_in, number_of_particles_in, size_t(store_output), 0, std::vector<const ProposalKernel *>(), double(number_of_particles_in), proposal_is_evaluated_in, smcfixed_flag_in, sequencer_limit_is_fixed_in, transform_proposed_particles, results_name_in)
  {
    this->chunk_size = grain_size_in;
    proposal_in->set_proposal_parameters(&this->algorithm_parameters);

    std::vector<size_t> indices;
    indices.reserve(likelihood_estimators_in.size());
    for (size_t i = 0; i < likelihood_estimators_in.size(); ++i)
      indices.push_back(i);
    this->index.reset(new VectorIndex(indices));

    this->factors.reset(new VectorFactors(likelihood_estimators_in));

    // Need to construct LikelihoodEstimator to read in to this constructor.
    this->particle_simulator.reset(new ParameterParticleSimulator(proposal_in,
                                                              likelihood_estimators_in));

    this->parallel = parallel_in;
    this->the_worker.reset(new SequentialSMCWorker(this));

    std::vector<double> schedule_in;
    schedule_in.push_back(target_value_in);
    SMCCriterion *smc_criterion = new PositiveSMCCriterion();
    this->sequencer = Sequencer(this->the_worker.get(),
                                schedule_in,
                                target_variable_in,
                                25,
                                smc_criterion);

    // sort out ABC scaling

    size_t data_size = 0;
    for (auto i = measurement_variables_in.begin();
         i != measurement_variables_in.end();
         ++i)
    {
      data_size = data_size + (*this->data)[*i].n_elem;
    }

    // if scale variable is not set, we are not using a scale at all, and we set it to be ones
    if (scale_variable_in == "")
    {
      this->sequencer.scale_variable = "default_scale";
      arma::colvec scale = arma::colvec(data_size);
      scale.fill(1.0);
      this->sequencer.schedule_parameters[this->sequencer.scale_variable] = scale;
      this->sequencer.find_scale = false;
    }
    else // if the scale variable is set, then either we have read in a valid scale, or it must be estimated through simulation
    {
      this->sequencer.scale_variable = scale_variable_in;

      if (scale_in.n_elem == data_size)
      {
        this->sequencer.schedule_parameters[this->sequencer.scale_variable] = scale_in;
        this->sequencer.find_scale = false;
      }
      else
      {
        Rcpp::stop("Scale variable is set, but the scale is not the correct size (it does not match the dimensions of the observed data).");
      }
    }
  }

  // Copy constructor for the ImportanceSampler class.
  ImportanceSampler::ImportanceSampler(const ImportanceSampler &another)
      : SMC(another)
  {
    this->make_copy(another);
  }

  // Destructor for the ImportanceSampler class.
  ImportanceSampler::~ImportanceSampler()
  {
    this->index.reset();
  }

  void ImportanceSampler::operator=(const ImportanceSampler &another)
  {
    if (this == &another)
    { // if a==a
      return;
    }

    this->index.reset();

    SMC::operator=(another);
    this->make_copy(another);
  }

  // Set the sequence variable (e.g. the ABC tolerance) to the value the
  // caller passes in parameters. The particles see it because they share the
  // sequencer's schedule_parameters (assignment writes into the shared matrix).
  void ImportanceSampler::set_target_from_parameters(const Parameters &parameters)
  {
    const std::string &variable = this->sequencer.variable_names.back();
    // (stops with a message naming the variable if the caller does not set it)
    double value = parameters[variable][0];
    this->sequencer.schedule_parameters[variable] = value;
  }

  void ImportanceSampler::set_abc_scale(Particles *particles)
  {
    if (this->sequencer.find_scale == true)
    {
      arma::mat packed_members;
      std::vector<arma::rowvec> partially_packed_members_row;
      // can be done in parallel
      for (auto i = particles->particles.begin();
           i != particles->particles.end();
           ++i)
      {
        arma::colvec packed_parameters = (*i)->back().get_colvec(this->sequencer.scale_variables);
        partially_packed_members_row.push_back(arma::conv_to<arma::rowvec>::from(packed_parameters));
      }

      // must be serial
      for (size_t i = 0;
           i < partially_packed_members_row.size();
           ++i)
      {
        if (i == 0)
        {
          packed_members = partially_packed_members_row[i];
        }
        else
        {
          packed_members = join_cols(packed_members, partially_packed_members_row[i]);
        }
      }

      arma::colvec abc_scale = arma::stddev(packed_members, 0, 0);
      for (size_t i = 0;
           i < abc_scale.n_elem;
           ++i)
      {
        if (abc_scale[i] == 0.0)
          abc_scale[i] = 1.0;
      }
      this->sequencer.schedule_parameters[this->sequencer.scale_variable] = abc_scale;
    }
  }

  SMC *ImportanceSampler::smc_duplicate() const
  {
    return (new ImportanceSampler(*this));
  }

  LikelihoodEstimator *ImportanceSampler::duplicate() const
  {
    return (new ImportanceSampler(*this));
  }

  void ImportanceSampler::make_copy(const ImportanceSampler &another)
  {
    this->follows_caller_target = another.follows_caller_target;
    if (another.index != NULL)
      this->index.reset(another.index->duplicate());
    else
      this->index.reset();
  }

  SMCOutput *ImportanceSampler::specific_run()
  {
    SMCOutput *simulation = this->initialise_smc();
    ProgressScope progress(this->show_progress, double(this->number_of_particles), "Importance sampling");
    if (this->parallel && !this->can_run_in_chunks())
      this->parallel_refused = true;
    if (this->can_run_in_chunks())
    {
      this->progress_counts_points = this->show_progress;
      this->simulate_and_weight_in_chunks(simulation);
      this->progress_counts_points = false;
      this->evaluate_smcadaptive_part_given_smcfixed_smc(simulation);
    }
    else
    {
      this->simulate_smc(simulation);
      this->evaluate_smc(simulation);
    }
    simulation->normalise_and_resample_weights();
    return simulation;
  }

  bool ImportanceSampler::can_run_in_chunks() const
  {
    // The weights computed while a chunk is simulated must be the final ones,
    // so the target must be fixed: one schedule point and a criterion that
    // always accepts it (no bisection), and no scale computed from all the
    // points (ABC). Weighting at other targets is G2 (adaptive summaries).
    if (this->force_rich_generations || this->proposed_particles_inputted || this->follows_caller_target)
      return false;
    if ((this->particle_simulator == NULL) || (this->number_of_particles == 0))
      return false;
    if (dynamic_cast<const SequentialSMCWorker*>(this->the_worker.get()) == NULL)
      return false;
    if (!this->sequencer_limit_is_fixed || this->sequencer.find_scale)
      return false;
    if ((this->sequencer.schedules.size() != 1) || (this->sequencer.schedules.back().size() != 1))
      return false;
    if (dynamic_cast<const PositiveSMCCriterion*>(this->sequencer.criterion.get()) == NULL)
      return false;
    return true;
  }

  void ImportanceSampler::simulate_and_weight_in_chunks(SMCOutput *current_state)
  {
    this->simulate_in_chunks(current_state, NULL, false);
  }

  bool ImportanceSampler::can_run_flat_conditioned() const
  {
    // Moving the smcfixed part to simulation time leaves every random number
    // unchanged only if that part uses none (checked per estimator), and
    // weighting from summaries needs every factor of the index to have one.
    if (this->force_rich_generations || this->proposed_particles_inputted)
      return false;
    if ((this->particle_simulator == NULL) || (this->number_of_particles == 0))
      return false;
    if (dynamic_cast<const SequentialSMCWorker*>(this->the_worker.get()) == NULL)
      return false;
    if (this->sequencer.find_scale)
      return false;
    if ((this->factors == NULL) || !this->factors->supports_flat_generations(this->index.get()))
      return false;
    return true;
  }

  bool ImportanceSampler::defers_simulation() const
  {
    return this->sequencer_limit_is_fixed && !this->proposal_is_evaluated && this->can_run_flat_conditioned();
  }

  // (dev/plan.md 3.21) the estimate is the mean of N weights, each at most
  // M = the sum of the factors' maxima (an ABC kernel's largest value): after
  // n points with weights summing to S it is at most (S + (N - n) e^M)/N, so
  // the run stops when the log of that is below threshold. The points are
  // made in chunks with their own random-number streams, so a run that does
  // not stop is the run without the check, except for the chunk size, which
  // is made small enough (16 chunks) for the check to be useful.
  bool ImportanceSampler::simulate_smc_with_threshold(SMCOutput* simulation,
                                                      const Parameters &conditioned_on_parameters,
                                                      double threshold)
  {
    const Parameters &schedule_parameters = this->sequencer.schedule_parameters;
    double M = (this->factors != NULL) ? this->factors->log_maximum(this->index.get(), schedule_parameters) : arma::datum::inf;
    size_t N = this->number_of_particles;
    if (!std::isfinite(M) || (N == 0) || !this->can_run_flat_conditioned())
    {
      this->simulate_smc(simulation, conditioned_on_parameters);
      return false;
    }
    double log_N = std::log(double(N));
    double S = -arma::datum::inf;
    size_t done_before = 0;
    std::vector<double> values;
    const Index* index = this->index.get();
    auto log_add = [](double a, double b) {
      if (a == -arma::datum::inf) return b;
      if (b == -arma::datum::inf) return a;
      double top = std::max(a, b);
      return top + std::log(std::exp(a - top) + std::exp(b - top));
    };
    this->chunk_stop = [&](Particles &generation, size_t done) -> bool
    {
      size_t count = done - done_before;
      if ((count > 0) && (generation.adaptive_summary_evaluator != NULL))
      {
        values.resize(count);
        size_t stride = generation.adaptive_summaries.n_rows;
        generation.adaptive_summary_evaluator->smcadaptive_from_summaries(index,
                                                                         generation.adaptive_summaries.memptr() + done_before*stride,
                                                                         stride, count, schedule_parameters, values.data());
        for (size_t j = 0; j < count; ++j)
          S = log_add(S, values[j]);
      }
      done_before = done;
      double rest = (done < N) ? std::log(double(N - done)) + M : -arma::datum::inf;
      return log_add(S, rest) - log_N < threshold;
    };
    size_t chunk_size = this->chunk_size;
    size_t small = (N + 15)/16;
    if ((this->chunk_size == 0) || (this->chunk_size > small))
      this->chunk_size = small;
    this->simulate_smc(simulation, conditioned_on_parameters);
    this->chunk_size = chunk_size;
    this->chunk_stop = nullptr;
    return this->chunk_stopped;
  }

  void ImportanceSampler::simulate_in_chunks(SMCOutput *current_state,
                                             const Parameters* conditioned_on_parameters,
                                             bool summarise)
  {
    this->simulate_generation_in_chunks(current_state,
                                        conditioned_on_parameters,
                                        summarise,
                                        this->index.get(),
                                        this->lag > 0);
  }

  bool ImportanceSampler::can_run_in_parallel() const
  {
    // Chunks may run at the same time if nothing they call shares mutable
    // state, or if each group can have its own copy of what does (G5).
    return (this->factors != NULL) &&
      (this->factors->supports_parallel_chunks() || this->factors->supports_parallel_copies());
  }

  bool ImportanceSampler::supports_parallel_copies() const
  {
    return this->can_run_flat_conditioned() && (this->lag == 0) &&
      (this->factors != NULL) && this->factors->supports_parallel_copies();
  }

  SMCOutput *ImportanceSampler::specific_initialise_smc()
  {
    SMCOutput *output = new SMCOutput(this, this->lag, this->lag_proposed, this->results_name);
    return output;
  }

  void ImportanceSampler::simulate_smc(SMCOutput *current_state)
  {
    this->simulate_proposal(current_state);
  }

  void ImportanceSampler::evaluate_smc(SMCOutput *current_state)
  {
    //

    this->evaluate_smcfixed_part_smc(current_state);
    this->evaluate_smcadaptive_part_given_smcfixed_smc(current_state);
  }

  void ImportanceSampler::evaluate_smcfixed_part_smc(SMCOutput *current_state)
  {

    this->the_worker->smcfixed_weight(this->index.get(),
                                      current_state->back());

  }

  void ImportanceSampler::evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput *current_state)
  {
    if (this->follows_caller_target)
    {
      this->weight_at_caller_target(current_state);
      return;
    }

    this->set_abc_scale(&current_state->back());

    this->sequencer.find_next_target_bisection(current_state, this->index.get());

    current_state->log_likelihood = current_state->calculate_latest_log_normalising_constant_ratio();
    current_state->set_llhd(current_state->log_likelihood);

    current_state->back().schedule_parameters = this->sequencer.schedule_parameters;

    current_state->set_time_and_reset_start();

    current_state->terminate();
  }

  MoveOutput *ImportanceSampler::move(RandomNumberGenerator &rng,
                                      const Particle &particle) const
  {
    return new SinglePointMoveOutput(std::move(particle));
  }

  void ImportanceSampler::weight_for_adapting_sequence(const Index *index,
                                                       Particles &current_particles)
  {
    if (current_particles.is_flat())
    {
      if (current_particles.adaptive_summary_state == 1)
      {
        // from the adaptive summaries kept when it was simulated
        this->weight_from_summaries_if_possible(index, current_particles);
      }
      else
      {
        // weights computed while the generation was simulated (at its only
        // target: see can_run_in_chunks)
        this->the_worker->set_unnormalised_log_incremental_weights(current_particles.flat_incremental_log_weights);
      }
      return;
    }
    this->the_worker->smcadaptive_given_smcfixed_weight(index,
                                                        current_particles);
  }

  // (dev/plan.md 3.21, step (b)) the points were drawn from the prior of the
  // variable integrated out, which is the proposal and cancels, so each
  // point's weight at the current target is its value of the other factors
  // (the data's densities, which read the caller's target, e.g. an ABC
  // tolerance, through the shared parameters): the estimate is their mean.
  // The points are kept (no terminate()), to be weighted again at the
  // caller's next target.
  void ImportanceSampler::weight_at_caller_target(SMCOutput *current_state)
  {
    Particles &generation = current_state->back();
    size_t n = generation.size();
    double log_likelihood = -arma::datum::inf;
    if (n > 0)
    {
      arma::colvec values(n);
      for (size_t i = 0; i < n; ++i)
        values[i] = generation[i]->back().evaluate_smcadaptive_part_given_smcfixed_likelihoods(this->index.get());
      log_likelihood = log_sum_exp(values) - std::log(double(n));
    }
    current_state->log_likelihood = log_likelihood;
    current_state->set_llhd(log_likelihood);
  }

  SMCOutput *ImportanceSampler::specific_run(const Parameters &conditioned_on_parameters)
  {
    SMCOutput *simulation = this->initialise_smc(conditioned_on_parameters);
    this->simulate_smc(simulation, conditioned_on_parameters);
    this->evaluate_smc(simulation, conditioned_on_parameters);
    simulation->normalise_and_resample_weights();
    return simulation;
  }

  SMCOutput *ImportanceSampler::specific_initialise_smc(const Parameters &conditioned_on_parameters)
  {
    SMCOutput *output = new SMCOutput(this, this->lag, this->lag_proposed, this->results_name);
    return output;
  }

  void ImportanceSampler::simulate_smc(SMCOutput *current_state,
                                       const Parameters &conditioned_on_parameters)
  {
    if (this->can_run_flat_conditioned())
    {
      this->simulate_in_chunks(current_state, &conditioned_on_parameters, true);
    }
    else
    {
      if (this->parallel)
        this->parallel_refused = true;
      this->simulate_proposal(current_state, conditioned_on_parameters);
    }
  }

  void ImportanceSampler::evaluate_smc(SMCOutput *current_state,
                                       const Parameters &conditioned_on_parameters)
  {

    this->evaluate_smcfixed_part_smc(current_state, conditioned_on_parameters);
    this->evaluate_smcadaptive_part_given_smcfixed_smc(current_state, conditioned_on_parameters);
  }

  void ImportanceSampler::evaluate_smcfixed_part_smc(SMCOutput *current_state,
                                                     const Parameters &conditioned_on_parameters)
  {
    // a flat generation's smcfixed parts were evaluated when it was simulated
    if (current_state->back().is_flat())
      return;

    this->the_worker->smcfixed_weight(this->index.get(),
                                      current_state->back());
  }

  void ImportanceSampler::evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput *current_state,
                                                                       const Parameters &conditioned_on_parameters)
  {
    if (this->follows_caller_target)
    {
      this->weight_at_caller_target(current_state);
      return;
    }
    if (!this->sequencer_limit_is_fixed)
    {
      // The target (e.g. the ABC tolerance) is set by the caller, through
      // conditioned_on_parameters. An importance sampler has a single target,
      // so it is reweighted directly at that value: no bisection, and no state
      // carried over from previous calls (this object is shared by every
      // particle of the calling algorithm).
      this->set_target_from_parameters(conditioned_on_parameters);
      this->set_abc_scale(&current_state->back());
      this->weight_for_adapting_sequence(this->index.get(),
                                         current_state->back());
      current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    }
    else
    {
      this->sequencer.reset();
      this->set_abc_scale(&current_state->back());
      this->sequencer.find_next_target_bisection(current_state, this->index.get());
    }

    current_state->log_likelihood = current_state->calculate_latest_log_normalising_constant_ratio();
    current_state->set_llhd(current_state->log_likelihood);

    current_state->back().schedule_parameters = this->sequencer.schedule_parameters;

    current_state->set_time_and_reset_start();

    // When the target is set by the caller (e.g. an ABC tolerance chosen by an
    // outer SMC), the simulations are kept: the caller reweights them at
    // new targets. terminate() would discard every stored particle set when
    // no output is stored.
    if (this->sequencer_limit_is_fixed)
      current_state->terminate();
  }

  void ImportanceSampler::subsample_simulate_smc(SMCOutput *current_state)
  {
    this->simulate_proposal(current_state);
  }

  void ImportanceSampler::subsample_simulate_smc(SMCOutput *current_state,
                                                 const Parameters &conditioned_on_parameters)
  {
    this->simulate_proposal(current_state, conditioned_on_parameters);
  }

  void ImportanceSampler::subsample_evaluate_smc(SMCOutput *current_state)
  {

    this->subsample_evaluate_smcfixed_part_smc(current_state);
    this->subsample_evaluate_smcadaptive_part_given_smcfixed_smc(current_state);
  }

  void ImportanceSampler::subsample_evaluate_smc(SMCOutput *current_state,
                                                 const Parameters &conditioned_on_parameters)
  {

    this->subsample_evaluate_smcfixed_part_smc(current_state, conditioned_on_parameters);
    this->subsample_evaluate_smcadaptive_part_given_smcfixed_smc(current_state, conditioned_on_parameters);
  }

  void ImportanceSampler::subsample_evaluate_smcfixed_part_smc(SMCOutput *current_state)
  {

    this->the_worker->subsample_smcfixed_weight(this->index.get(),
                                                current_state->back());
  }

  void ImportanceSampler::subsample_evaluate_smcfixed_part_smc(SMCOutput *current_state,
                                                               const Parameters &conditioned_on_parameters)
  {

    this->the_worker->subsample_smcfixed_weight(this->index.get(),
                                                current_state->back());
  }

  void ImportanceSampler::subsample_evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput *current_state)
  {
    // set sequencer to have values from conditioned_on_parameters

    this->set_abc_scale(&current_state->back());

    this->sequencer.subsample_find_next_target_bisection(current_state, this->index.get());

    // A subsample run returns subsample_log_likelihood (see LikelihoodEstimatorOutput::subsample_evaluate).
    current_state->subsample_log_likelihood = current_state->calculate_latest_log_normalising_constant_ratio();

    current_state->back().schedule_parameters = this->sequencer.schedule_parameters;

    current_state->set_llhd(current_state->subsample_log_likelihood);
    current_state->set_time_and_reset_start();

    current_state->terminate();
  }

  void ImportanceSampler::subsample_evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput *current_state,
                                                                                 const Parameters &conditioned_on_parameters)
  {
    if (!this->sequencer_limit_is_fixed)
    {
      // As in evaluate_smcadaptive_part_given_smcfixed_smc: single target set
      // by the caller.
      this->set_target_from_parameters(conditioned_on_parameters);
      this->set_abc_scale(&current_state->back());
      this->the_worker->subsample_smcadaptive_given_smcfixed_weight(this->index.get(),
                                                                    current_state->back());
      current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    }
    else
    {
      this->sequencer.reset();
      this->set_abc_scale(&current_state->back());
      this->sequencer.subsample_find_next_target_bisection(current_state, this->index.get());
    }

    // A subsample run returns subsample_log_likelihood (see LikelihoodEstimatorOutput::subsample_evaluate).
    current_state->subsample_log_likelihood = current_state->calculate_latest_log_normalising_constant_ratio();

    current_state->back().schedule_parameters = this->sequencer.schedule_parameters;

    current_state->set_llhd(current_state->subsample_log_likelihood);
    current_state->set_time_and_reset_start();

    // When the target is set by the caller (e.g. an ABC tolerance chosen by an
    // outer SMC), the simulations are kept: the caller reweights them at
    // new targets. terminate() would discard every stored particle set when
    // no output is stored.
    if (this->sequencer_limit_is_fixed)
      current_state->terminate();
  }

  MoveOutput *ImportanceSampler::subsample_move(RandomNumberGenerator &rng,
                                                const Particle &particle) const
  {
    return new SinglePointMoveOutput(std::move(particle));
  }

  void ImportanceSampler::subsample_weight_for_adapting_sequence(const Index *index,
                                                                 Particles &current_particles)
  {
    this->the_worker->subsample_smcadaptive_given_smcfixed_weight(index,
                                                                  current_particles);
  }

  // Comment for later...
  // IS: simulate
  // SMC w MCMC: t=1 simulate, t>1 loop resample-move, then weight, until stopping point reached with resample-move being the final step
  // PF/SMC: t=1 simulate, t>1 resample, sim prop, then weight, until stopping point reached with sim prop being final step
  // void ImportanceSampler::smc_simulate(SMCOutput* current_state)
  //{
  // this->output->add_proposed_particles(particles);

}
