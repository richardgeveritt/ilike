#include "importance_sampler.h"
#include "smc_worker.h"
// #include "rcppparallel_smc_worker.h"
#include "sequential_smc_worker.h"
#include "smc_output.h"
#include "exact_likelihood_estimator.h"
#include "parameter_particle_simulator.h"
#include "move_output.h"
#include "single_point_move_output.h"
#include "independent_proposal_kernel.h"
// #include "custom_independent_proposal_kernel.h"
#include "custom_distribution_proposal_kernel.h"
#include "vector_factors.h"
#include "vector_index.h"
#include "positive_smc_criterion.h"
#include "factor_variables.h"
#include "particle_simulator.h"
#include <memory>
// #include "rcppparallel_smc_worker.h"

namespace ilike
{
  ImportanceSampler::ImportanceSampler()
      : SMC()
  {
    this->index = NULL;
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
    // Parameters candidate_parameters = proposal->independent_simulate(*this->rng);

    std::vector<size_t> indices;
    std::vector<LikelihoodEstimator *> likelihood_estimators;
    likelihood_estimators.reserve(1);
    likelihood_estimators.push_back(new ExactLikelihoodEstimator(rng_in,
                                                                 seed_in,
                                                                 data_in,
                                                                 evaluate_log_likelihood_in,
                                                                 true));
    indices.push_back(0);
    this->index = new VectorIndex(indices);

    this->factors = new VectorFactors(likelihood_estimators);

    // Need to construct LikelihoodEstimator to read in to this constructor.
    this->particle_simulator = new ParameterParticleSimulator(proposal,
                                                              likelihood_estimators);

    if (parallel_in == TRUE)
    {
      this->the_worker = NULL;
    }
    else
    {
      this->the_worker = new SequentialSMCWorker(this);
    }

    std::vector<double> schedule_in;
    schedule_in.push_back(0.0);
    schedule_in.push_back(1.0);
    SMCCriterion *smc_criterion = new PositiveSMCCriterion();
    this->sequencer = Sequencer(this->the_worker,
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
    // Parameters candidate_parameters = proposal->independent_simulate(*this->rng);

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
    this->index = new VectorIndex(indices);

    this->factors = new VectorFactors(likelihood_estimators);

    // Need to construct LikelihoodEstimator to read in to this constructor.
    this->particle_simulator = new ParameterParticleSimulator(proposal,
                                                              likelihood_estimators);

    /*
     for (auto i=likelihood_estimators.begin();
     i!=likelihood_estimators.end();
     ++i)
     {
     (*i)->setup(candidate_parameters);
     }
     */

    // this->model_and_algorithm.particle_simulator = new ParameterParticleSimulator(simulate_proposal_in,
    //                                                                               this->model_and_algorithm.likelihood_estimators,
    //                                                                               "u");

    if (parallel_in == TRUE)
    {
    }
    else
    {
      this->the_worker = new SequentialSMCWorker(this);
    }

    std::vector<double> schedule_in;
    schedule_in.push_back(0.0);
    schedule_in.push_back(1.0);
    SMCCriterion *smc_criterion = new PositiveSMCCriterion();
    this->sequencer = Sequencer(this->the_worker,
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
    // Parameters candidate_parameters = proposal_in->independent_simulate(*this->rng);
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
    this->index = new VectorIndex(indices);

    this->factors = new VectorFactors(likelihood_estimators);

    // Need to construct LikelihoodEstimator to read in to this constructor.
    this->particle_simulator = new ParameterParticleSimulator(proposal_in,
                                                              likelihood_estimators);

    // this->model_and_algorithm.particle_simulator = new ParameterParticleSimulator(simulate_proposal_in,
    //                                                                               this->model_and_algorithm.likelihood_estimators,
    //                                                                               "u");

    if (parallel_in == TRUE)
    {
    }
    else
    {
      this->the_worker = new SequentialSMCWorker(this);
    }

    std::vector<double> schedule_in;
    schedule_in.push_back(0.0);
    schedule_in.push_back(1.0);
    SMCCriterion *smc_criterion = new PositiveSMCCriterion();
    this->sequencer = Sequencer(this->the_worker,
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

    // Parameters candidate_parameters = proposal_in->independent_simulate(*this->rng);
    //std::vector<size_t> indices;
    //indices.reserve(likelihood_estimators_in.size());
    //for (size_t i = 0; i < likelihood_estimators_in.size(); ++i)
    //  indices.push_back(i);
    //this->index = new VectorIndex(indices);
    this->index = without_cancelled_index;

    this->factors = new VectorFactors(likelihood_estimators_in);

    /*
     for (auto i=likelihood_estimators_in.begin();
     i!=likelihood_estimators_in.end();
     ++i)
     {
     (*i)->setup(candidate_parameters);
     }
     */

    // Need to construct LikelihoodEstimator to read in to this constructor.
    this->particle_simulator = new ParameterParticleSimulator(proposal_in,
                                                              likelihood_estimators_in);

    // this->model_and_algorithm.particle_simulator = new ParameterParticleSimulator(simulate_proposal_in,
    //                                                                               this->model_and_algorithm.likelihood_estimators,
    //                                                                               "u");

    if (parallel_in == TRUE)
    {
      // this->the_worker = new RcppParallelSMCWorker(this,grain_size_in);
      this->the_worker = NULL;
    }
    else
    {
      this->the_worker = new SequentialSMCWorker(this);
    }

    std::vector<double> schedule_in;
    schedule_in.push_back(1.0);
    SMCCriterion *smc_criterion = new PositiveSMCCriterion();
    this->sequencer = Sequencer(this->the_worker,
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

    // Parameters candidate_parameters = proposal_in->independent_simulate(*this->rng);
    //std::vector<size_t> indices;
    //indices.reserve(likelihood_estimators_in.size());
    //for (size_t i = 0; i < likelihood_estimators_in.size(); ++i)
    //  indices.push_back(i);
    //this->index = new VectorIndex(indices);
    this->index = without_cancelled_index;

    this->factors = new VectorFactors(likelihood_estimators_in);

    /*
     for (auto i=likelihood_estimators_in.begin();
     i!=likelihood_estimators_in.end();
     ++i)
     {
     (*i)->setup(candidate_parameters);
     }
     */

    // Need to construct LikelihoodEstimator to read in to this constructor.
    this->particle_simulator = new ParameterParticleSimulator(proposal_in,
                                                              likelihood_estimators_in);

    // this->model_and_algorithm.particle_simulator = new ParameterParticleSimulator(simulate_proposal_in,
    //                                                                               this->model_and_algorithm.likelihood_estimators,
    //                                                                               "u");

    if (parallel_in == TRUE)
    {
      // this->the_worker = new RcppParallelSMCWorker(this,grain_size_in);
      this->the_worker = NULL;
    }
    else
    {
      this->the_worker = new SequentialSMCWorker(this);
    }

    std::vector<double> schedule_in;
    schedule_in.push_back(target_value_in);
    SMCCriterion *smc_criterion = new PositiveSMCCriterion();
    this->sequencer = Sequencer(this->the_worker,
                                schedule_in,
                                target_variable_in,
                                25,
                                smc_criterion);
    // this->sequencer_parameters = &this->sequencer.schedule_parameters;
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

    // Parameters candidate_parameters = proposal_in->independent_simulate(*this->rng);
    std::vector<size_t> indices;
    indices.reserve(likelihood_estimators_in.size());
    for (size_t i = 0; i < likelihood_estimators_in.size(); ++i)
      indices.push_back(i);
    this->index = new VectorIndex(indices);

    this->factors = new VectorFactors(likelihood_estimators_in);

    /*
     for (auto i=likelihood_estimators_in.begin();
     i!=likelihood_estimators_in.end();
     ++i)
     {
     (*i)->setup(candidate_parameters);
     }
     */

    // Need to construct LikelihoodEstimator to read in to this constructor.
    this->particle_simulator = new ParameterParticleSimulator(proposal_in,
                                                              likelihood_estimators_in);

    // this->model_and_algorithm.particle_simulator = new ParameterParticleSimulator(simulate_proposal_in,
    //                                                                               this->model_and_algorithm.likelihood_estimators,
    //                                                                               "u");

    if (parallel_in == TRUE)
    {
      // this->the_worker = new RcppParallelSMCWorker(this,grain_size_in);
      this->the_worker = NULL;
    }
    else
    {
      this->the_worker = new SequentialSMCWorker(this);
    }

    std::vector<double> schedule_in;
    schedule_in.push_back(1.0);
    SMCCriterion *smc_criterion = new PositiveSMCCriterion();
    this->sequencer = Sequencer(this->the_worker,
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
        // this->sequencer.scale_variables = measurement_variables_in;
        // this->sequencer.find_scale = true;
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

    // Parameters candidate_parameters = proposal_in->independent_simulate(*this->rng);
    std::vector<size_t> indices;
    indices.reserve(likelihood_estimators_in.size());
    for (size_t i = 0; i < likelihood_estimators_in.size(); ++i)
      indices.push_back(i);
    this->index = new VectorIndex(indices);

    this->factors = new VectorFactors(likelihood_estimators_in);

    /*
     for (auto i=likelihood_estimators_in.begin();
     i!=likelihood_estimators_in.end();
     ++i)
     {
     (*i)->setup(candidate_parameters);
     }
     */

    // Need to construct LikelihoodEstimator to read in to this constructor.
    this->particle_simulator = new ParameterParticleSimulator(proposal_in,
                                                              likelihood_estimators_in);

    // this->model_and_algorithm.particle_simulator = new ParameterParticleSimulator(simulate_proposal_in,
    //                                                                               this->model_and_algorithm.likelihood_estimators,
    //                                                                               "u");

    if (parallel_in == TRUE)
    {
      // this->the_worker = new RcppParallelSMCWorker(this,grain_size_in);
      this->the_worker = NULL;
    }
    else
    {
      this->the_worker = new SequentialSMCWorker(this);
    }

    std::vector<double> schedule_in;
    schedule_in.push_back(target_value_in);
    SMCCriterion *smc_criterion = new PositiveSMCCriterion();
    this->sequencer = Sequencer(this->the_worker,
                                schedule_in,
                                target_variable_in,
                                25,
                                smc_criterion);
    // this->sequencer_parameters = &this->sequencer.schedule_parameters;

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
        // this->sequencer.scale_variables = measurement_variables_in;
        // this->sequencer.find_scale = true;
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
    if (this->index != NULL)
      delete this->index;
  }

  void ImportanceSampler::operator=(const ImportanceSampler &another)
  {
    if (this == &another)
    { // if a==a
      return;
    }

    if (this->index != NULL)
      delete this->index;

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
    if (another.index != NULL)
      this->index = another.index->duplicate();
    else
      this->index = NULL;
  }

  // void ImportanceSampler::smc_step()
  // {
  // }
  //
  // void ImportanceSampler::weight_update()
  // {
  // }

  SMCOutput *ImportanceSampler::specific_run()
  {
    SMCOutput *simulation = this->initialise_smc();
    if (this->can_run_in_chunks())
    {
      this->simulate_and_weight_in_chunks(simulation);
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
    if (this->force_rich_generations || this->proposed_particles_inputted)
      return false;
    if ((this->particle_simulator == NULL) || (this->number_of_particles == 0))
      return false;
    if (dynamic_cast<const SequentialSMCWorker*>(this->the_worker) == NULL)
      return false;
    if (!this->sequencer_limit_is_fixed || this->sequencer.find_scale)
      return false;
    if ((this->sequencer.schedules.size() != 1) || (this->sequencer.schedules.back().size() != 1))
      return false;
    if (dynamic_cast<const PositiveSMCCriterion*>(this->sequencer.criterion) == NULL)
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
    if (dynamic_cast<const SequentialSMCWorker*>(this->the_worker) == NULL)
      return false;
    if (this->sequencer.find_scale)
      return false;
    if ((this->factors == NULL) || !this->factors->supports_flat_generations(this->index))
      return false;
    return true;
  }

  void ImportanceSampler::simulate_in_chunks(SMCOutput *current_state,
                                             const Parameters* conditioned_on_parameters,
                                             bool summarise)
  {
    // What simulate_smc + evaluate_smcfixed_part_smc + the weighting step of
    // the bisection do (SequentialSMCWorker::specific_simulate,
    // SMCWorker::simulate, smcfixed_weight, smcadaptive_given_smcfixed_weight),
    // B points at a time. Per chunk the phases keep their order (simulate all
    // B, then the smcfixed parts, then the weights), so with one chunk
    // (chunk_size >= N) every call happens in the same order as before and
    // the results are identical. Chunk k uses its own random-number stream:
    // the stream of chunk 0 (the generator seeded as before) jumped k times
    // (2^192 steps each), so streams never overlap and depend only on the
    // seed and k.
    SMCWorker* worker = this->the_worker;
    size_t N = this->number_of_particles;
    size_t B = ((this->chunk_size == 0) || (this->chunk_size > N)) ? N : this->chunk_size;

    Particles* generation = current_state->add_particles();
    generation->flat = true;
    generation->flat_size = N;
    if (!summarise)
      generation->flat_incremental_log_weights.set_size(N);

    // What is kept of each point: its variables and factor outputs only if
    // the output stores generations (lag 0: an estimator inside a factor,
    // whose generations are never written).
    bool keep_points = (this->lag > 0);
    bool keep_factor_output = keep_points && this->write_factor_output;
    int summary_size = -1;
    bool compact_records = false;
    std::vector<double> records;

    std::vector<Particle> slots(B);
    RandomNumberGenerator chunk_start(*worker->get_rng());
    chunk_start.seed(worker->get_seed());

    size_t point = 0;
    for (size_t first = 0, k = 0; first < N; first += B, ++k)
    {
      size_t n = std::min(B, N - first);
      if (k > 0)
      {
        std::unique_ptr<dqrng::random_64bit_generator> next = chunk_start.clone(1);
        chunk_start = *dynamic_cast<RandomNumberGenerator*>(next.get());
      }
      RandomNumberGenerator local_rng(chunk_start);

      for (size_t j = 0; j < n; ++j)
      {
        if (conditioned_on_parameters == NULL)
          this->particle_simulator->simulate_into(local_rng,
                                                  slots[j],
                                                  this->factors,
                                                  &this->proposals_to_transform_for,
                                                  &this->proposals_to_find_gradient_for,
                                                  this->sequencer.schedule_parameters);
        else
          this->particle_simulator->simulate_into(local_rng,
                                                  slots[j],
                                                  this->factors,
                                                  &this->proposals_to_transform_for,
                                                  &this->proposals_to_find_gradient_for,
                                                  *conditioned_on_parameters,
                                                  this->sequencer.schedule_parameters);
        if (this->proposal_is_evaluated == true)
          slots[j].previous_target_evaluated = this->particle_simulator->evaluate(slots[j]);
      }

      if (k == 0)
      {
        // the seed bookkeeping and resampling variables of SMCWorker::simulate,
        // at the point where it happened before when there is one chunk
        worker->set_seed(worker->get_seed() + N);
        this->rng->seed(worker->get_seed());
        generation->simulate_resampling_variables(*this->rng);
        worker->set_seed(worker->get_seed() + N);
      }

      for (size_t j = 0; j < n; ++j)
        slots[j].evaluate_smcfixed_part_of_likelihoods(this->index);

      if (summarise)
      {
        for (size_t j = 0; j < n; ++j)
        {
          const FactorVariables* factor_variables = slots[j].factor_variables;
          int size = (factor_variables == NULL) ? -1 : factor_variables->adaptive_summary_size(this->index);
          if (first + j == 0)
          {
            if (size < 0)
              Rcpp::stop("ImportanceSampler::simulate_in_chunks: a factor has no adaptive summary.");
            summary_size = size;
            generation->adaptive_summaries.set_size(size_t(summary_size), N);
            generation->adaptive_previous_targets.set_size(N);
            generation->flat_summary_evaluator = factor_variables->duplicate();
            generation->flat_summary_evaluator->set_particle(NULL);
            generation->adaptive_summary_evaluator = generation->flat_summary_evaluator;
            generation->adaptive_summary_index = this->index;
            generation->adaptive_summary_state = 1;
          }
          else if (size != summary_size)
          {
            Rcpp::stop("ImportanceSampler::simulate_in_chunks: the points' adaptive summaries differ in size.");
          }
          if (summary_size > 0)
            factor_variables->adaptive_summary(this->index, generation->adaptive_summaries.colptr(first + j));
          generation->adaptive_previous_targets[first + j] = slots[j].previous_target_evaluated;
        }
      }
      else
      {
        for (size_t j = 0; j < n; ++j)
        {
          double a = slots[j].evaluate_smcadaptive_part_given_smcfixed_likelihoods(this->index);
          double b = slots[j].previous_target_evaluated;
          generation->flat_incremental_log_weights[first + j] = (a == -arma::datum::inf) ? -arma::datum::inf : a - b;
        }
      }

      for (size_t j = 0; (j < n) && keep_points; ++j, ++point)
      {
        Particle &slot = slots[j];
        if (point == 0)
        {
          if (keep_factor_output && (slot.factor_variables != NULL))
          {
            compact_records = slot.factor_variables->has_compact_values();
            if (compact_records)
            {
              generation->flat_points.set_number_of_extra_columns(slot.factor_variables->number_of_compact_values());
              records.resize(slot.factor_variables->number_of_compact_values());
              generation->flat_record_writer = slot.factor_variables->duplicate();
              generation->flat_record_writer->set_particle(NULL);
            }
            else
            {
              generation->flat_factor_variables.reserve(N);
            }
          }
          generation->flat_points.reserve_after_first(N);
        }
        if (compact_records)
        {
          slot.factor_variables->compact_values(records.data());
          generation->flat_points.append(slot.parameters, records.data());
        }
        else
        {
          generation->flat_points.append(slot.parameters);
          if (keep_factor_output && (slot.factor_variables != NULL))
          {
            // keep this point's factor outputs whole; the slot gets new ones
            slot.factor_variables->set_particle(NULL);
            generation->flat_factor_variables.push_back(slot.factor_variables);
            slot.factor_variables = NULL;
          }
        }
      }
    }

    generation->initialise_weights();
  }

  /*
   SMCOutput* ImportanceSampler::specific_run(const std::string &directory_name)
   {
   SMCOutput* simulation = this->initialise_smc();
   this->simulate_smc(simulation);
   this->evaluate_smc(simulation);
   simulation->normalise_and_resample_weights();
   simulation->write(directory_name);
   return simulation;
   }
   */

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
    /*
     //this->weight(current_state, conditioned_on_parameters);
     this->the_worker->weight(this->index,
     current_state->back());
     //current_state->initialise_next_step();
     current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());

     //this->log_likelihood = current_state->log_likelihood_pre_last_step + current_state->latest_log_normalising_constant_ratio();

     current_state->log_likelihood = current_state->calculate_latest_log_normalising_constant_ratio();

     if (this->sequencer_parameters!=NULL)
     current_state->back().schedule_parameters = *this->sequencer_parameters;
     */

    this->evaluate_smcfixed_part_smc(current_state);
    this->evaluate_smcadaptive_part_given_smcfixed_smc(current_state);
  }

  void ImportanceSampler::evaluate_smcfixed_part_smc(SMCOutput *current_state)
  {
    // this->smcfixed_weight(current_state, conditioned_on_parameters);

    this->the_worker->smcfixed_weight(this->index,
                                      current_state->back());

    /*
     if (this->sequencer_parameters!=NULL)
     {
     this->the_worker->smcfixed_weight(this->index,
     current_state->back(),
     *this->sequencer_parameters);
     }
     else
     {
     this->the_worker->smcfixed_weight(this->index,
     current_state->back());
     }
     */
    // current_state->initialise_next_step();
  }

  void ImportanceSampler::evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput *current_state)
  {
    // this->smcadaptive_given_smcfixed_weight(current_state, conditioned_on_parameters);
    // this->the_worker->smcadaptive_given_smcfixed_weight(this->index,
    //                                                     current_state->back());
    // current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());

    this->set_abc_scale(&current_state->back());

    this->sequencer.find_next_target_bisection(current_state, this->index);
    // this->log_likelihood = current_state->log_likelihood_pre_last_step + current_state->latest_log_normalising_constant_ratio();

    current_state->log_likelihood = current_state->calculate_latest_log_normalising_constant_ratio();
    current_state->set_llhd(current_state->log_likelihood);

    // if (this->sequencer_parameters!=NULL)
    //   current_state->back().schedule_parameters = *this->sequencer_parameters;

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

  /*
   void ImportanceSampler::weight_for_adapting_sequence(const Index* index,
   Particles &current_particles,
   const Parameters &conditioned_on_parameters)
   {
   this->the_worker->smcadaptive_given_smcfixed_weight(index,
   current_particles,
   conditioned_on_parameters);
   }
   */

  SMCOutput *ImportanceSampler::specific_run(const Parameters &conditioned_on_parameters)
  {
    SMCOutput *simulation = this->initialise_smc(conditioned_on_parameters);
    this->simulate_smc(simulation, conditioned_on_parameters);
    this->evaluate_smc(simulation, conditioned_on_parameters);
    simulation->normalise_and_resample_weights();
    return simulation;
  }

  /*
   SMCOutput* ImportanceSampler::specific_run(const std::string &directory_name,
   const Parameters &conditioned_on_parameters)
   {
   SMCOutput* simulation = this->initialise_smc(conditioned_on_parameters);
   this->simulate_smc(simulation, conditioned_on_parameters);
   this->evaluate_smc(simulation, conditioned_on_parameters);
   simulation->normalise_and_resample_weights();
   simulation->write(directory_name);
   return simulation;
   }
   */

  SMCOutput *ImportanceSampler::specific_initialise_smc(const Parameters &conditioned_on_parameters)
  {
    SMCOutput *output = new SMCOutput(this, this->lag, this->lag_proposed, this->results_name);
    return output;
  }

  void ImportanceSampler::simulate_smc(SMCOutput *current_state,
                                       const Parameters &conditioned_on_parameters)
  {
    if (this->can_run_flat_conditioned())
      this->simulate_in_chunks(current_state, &conditioned_on_parameters, true);
    else
      this->simulate_proposal(current_state, conditioned_on_parameters);
  }

  void ImportanceSampler::evaluate_smc(SMCOutput *current_state,
                                       const Parameters &conditioned_on_parameters)
  {
    /*
     //this->weight(current_state, conditioned_on_parameters);
     this->the_worker->weight(this->index,
     current_state->back(),
     conditioned_on_parameters);
     //current_state->initialise_next_step();
     current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());

     //this->log_likelihood = current_state->log_likelihood_pre_last_step + current_state->latest_log_normalising_constant_ratio();

     current_state->log_likelihood = current_state->calculate_latest_log_normalising_constant_ratio();

     if (this->sequencer_parameters!=NULL)
     current_state->back().schedule_parameters = *this->sequencer_parameters;
     */

    this->evaluate_smcfixed_part_smc(current_state, conditioned_on_parameters);
    this->evaluate_smcadaptive_part_given_smcfixed_smc(current_state, conditioned_on_parameters);
  }

  void ImportanceSampler::evaluate_smcfixed_part_smc(SMCOutput *current_state,
                                                     const Parameters &conditioned_on_parameters)
  {
    // a flat generation's smcfixed parts were evaluated when it was simulated
    if (current_state->back().is_flat())
      return;
    // this->smcfixed_weight(current_state, conditioned_on_parameters);
    /*
     if (this->sequencer_parameters!=NULL)
     {
     Parameters all_parameters = conditioned_on_parameters.merge(*this->sequencer_parameters);
     this->the_worker->smcfixed_weight(this->index,
     current_state->back(),
     all_parameters);
     }
     else
     {
     this->the_worker->smcfixed_weight(this->index,
     current_state->back(),
     conditioned_on_parameters);
     }
     */

    this->the_worker->smcfixed_weight(this->index,
                                      current_state->back());
  }

  void ImportanceSampler::evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput *current_state,
                                                                       const Parameters &conditioned_on_parameters)
  {
    if (!this->sequencer_limit_is_fixed)
    {
      // The target (e.g. the ABC tolerance) is set by the caller, through
      // conditioned_on_parameters. An importance sampler has a single target,
      // so it is reweighted directly at that value: no bisection, and no state
      // carried over from previous calls (this object is shared by every
      // particle of the calling algorithm).
      this->set_target_from_parameters(conditioned_on_parameters);
      this->set_abc_scale(&current_state->back());
      this->weight_for_adapting_sequence(this->index,
                                         current_state->back());
      current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    }
    else
    {
      this->sequencer.reset();
      this->set_abc_scale(&current_state->back());
      this->sequencer.find_next_target_bisection(current_state, this->index);
    }
    // this->the_worker->smcadaptive_given_smcfixed_weight(this->index,
    //                                                     current_state->back(),
    //                                                     conditioned_on_parameters);
    // current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());

    // this->log_likelihood = current_state->log_likelihood_pre_last_step + current_state->latest_log_normalising_constant_ratio();

    current_state->log_likelihood = current_state->calculate_latest_log_normalising_constant_ratio();
    current_state->set_llhd(current_state->log_likelihood);

    // if (this->sequencer_parameters!=NULL)
    //   current_state->back().schedule_parameters = *this->sequencer_parameters;
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

    /*
     this->the_worker->subsample_weight(this->index,
     current_state->back());
     current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
     current_state->subsample_log_likelihood = current_state->calculate_latest_log_normalising_constant_ratio();
     */

    this->subsample_evaluate_smcfixed_part_smc(current_state);
    this->subsample_evaluate_smcadaptive_part_given_smcfixed_smc(current_state);
  }

  void ImportanceSampler::subsample_evaluate_smc(SMCOutput *current_state,
                                                 const Parameters &conditioned_on_parameters)
  {

    /*
     this->the_worker->subsample_weight(this->index,
     current_state->back());
     current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());

     current_state->subsample_log_likelihood = current_state->calculate_latest_log_normalising_constant_ratio();
     */

    this->subsample_evaluate_smcfixed_part_smc(current_state, conditioned_on_parameters);
    this->subsample_evaluate_smcadaptive_part_given_smcfixed_smc(current_state, conditioned_on_parameters);
  }

  void ImportanceSampler::subsample_evaluate_smcfixed_part_smc(SMCOutput *current_state)
  {
    // this->smcfixed_weight(current_state, conditioned_on_parameters);
    /*
     if (this->sequencer_parameters!=NULL)
     {
     Parameters all_parameters = conditioned_on_parameters.merge(*this->sequencer_parameters);
     this->the_worker->subsample_smcfixed_weight(this->index,
     current_state->back(),
     all_parameters);
     }
     else
     {
     this->the_worker->subsample_smcfixed_weight(this->index,
     current_state->back(),
     conditioned_on_parameters);
     }
     */

    this->the_worker->subsample_smcfixed_weight(this->index,
                                                current_state->back());
    // current_state->initialise_next_step();
  }

  void ImportanceSampler::subsample_evaluate_smcfixed_part_smc(SMCOutput *current_state,
                                                               const Parameters &conditioned_on_parameters)
  {
    // this->smcfixed_weight(current_state, conditioned_on_parameters);
    /*
     if (this->sequencer_parameters!=NULL)
     {
     Parameters all_parameters = conditioned_on_parameters.merge(*this->sequencer_parameters);
     this->the_worker->subsample_smcfixed_weight(this->index,
     current_state->back(),
     all_parameters);
     }
     else
     {
     this->the_worker->subsample_smcfixed_weight(this->index,
     current_state->back(),
     conditioned_on_parameters);
     }
     */

    this->the_worker->subsample_smcfixed_weight(this->index,
                                                current_state->back());
    // current_state->initialise_next_step();
  }

  void ImportanceSampler::subsample_evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput *current_state)
  {
    // set sequencer to have values from conditioned_on_parameters

    this->set_abc_scale(&current_state->back());

    // this->sequencer.subsample_find_next_target_bisection(current_state,this->index,conditioned_on_parameters);
    this->sequencer.subsample_find_next_target_bisection(current_state, this->index);

    // this->smcadaptive_given_smcfixed_weight(current_state, conditioned_on_parameters);

    /*
     this->the_worker->subsample_smcadaptive_given_smcfixed_weight(this->index,
     current_state->back(),
     conditioned_on_parameters);
     current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
     */

    // this->log_likelihood = current_state->log_likelihood_pre_last_step + current_state->latest_log_normalising_constant_ratio();

    // A subsample run returns subsample_log_likelihood (see LikelihoodEstimatorOutput::subsample_evaluate).
    current_state->subsample_log_likelihood = current_state->calculate_latest_log_normalising_constant_ratio();

    // if (this->sequencer_parameters!=NULL)
    //   current_state->back().schedule_parameters = *this->sequencer_parameters;
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
      this->the_worker->subsample_smcadaptive_given_smcfixed_weight(this->index,
                                                                    current_state->back());
      current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    }
    else
    {
      this->sequencer.reset();
      this->set_abc_scale(&current_state->back());
      this->sequencer.subsample_find_next_target_bisection(current_state, this->index);
    }

    // this->smcadaptive_given_smcfixed_weight(current_state, conditioned_on_parameters);

    /*
     this->the_worker->subsample_smcadaptive_given_smcfixed_weight(this->index,
     current_state->back(),
     conditioned_on_parameters);
     current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
     */

    // this->log_likelihood = current_state->log_likelihood_pre_last_step + current_state->latest_log_normalising_constant_ratio();

    // A subsample run returns subsample_log_likelihood (see LikelihoodEstimatorOutput::subsample_evaluate).
    current_state->subsample_log_likelihood = current_state->calculate_latest_log_normalising_constant_ratio();

    // if (this->sequencer_parameters!=NULL)
    //   current_state->back().schedule_parameters = *this->sequencer_parameters;
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

  /*
   MoveOutput* ImportanceSampler::move(RandomNumberGenerator &rng,
   Particle &particle,
   const Parameters &conditioned_on_parameters)
   {
   return new SinglePointMoveOutput(std::move(particle));
   }
   */

  // void ImportanceSampler::weight_for_adapting_sequence(Particles &current_particles,
  //                                                      const Parameters &conditioned_on_parameters)
  //{
  //   this->the_worker->smcadaptive_given_smcfixed_weight(this->index,
  //                                                       current_particles,
  //                                                       conditioned_on_parameters);
  // }

  MoveOutput *ImportanceSampler::subsample_move(RandomNumberGenerator &rng,
                                                const Particle &particle) const
  {
    return new SinglePointMoveOutput(std::move(particle));
  }

  /*
   MoveOutput* ImportanceSampler::subsample_move(RandomNumberGenerator &rng,
   Particle &particle,
   const Parameters &conditioned_on_parameters)
   {
   return new SinglePointMoveOutput(std::move(particle));
   }
   */

  void ImportanceSampler::subsample_weight_for_adapting_sequence(const Index *index,
                                                                 Particles &current_particles)
  {
    this->the_worker->subsample_smcadaptive_given_smcfixed_weight(index,
                                                                  current_particles);
  }

  /*
   void ImportanceSampler::subsample_weight_for_adapting_sequence(const Index* index,
   Particles &current_particles,
   const Parameters &conditioned_on_parameters)
   {
   this->the_worker->subsample_smcadaptive_given_smcfixed_weight(index,
   current_particles,
   conditioned_on_parameters);
   }
   */

  /*
   void ImportanceSampler::weight(SMCOutput* current_state,
   const Parameters &conditioned_on_parameters)
   {
   this->the_worker->weight(conditioned_on_parameters);
   current_state->initialise_next_step();
   current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
   }

   void ImportanceSampler::smcfixed_weight(SMCOutput* current_state,
   const Parameters &conditioned_on_parameters)
   {
   this->the_worker->smcfixed_weight(conditioned_on_parameters);
   current_state->initialise_next_step();
   }


   void ImportanceSampler::smcadaptive_given_smcfixed_weight(SMCOutput* current_state,
   const Parameters &conditioned_on_parameters)
   {
   this->the_worker->smcadaptive_given_smcfixed_weight(conditioned_on_parameters);
   current_state->weight_update(this->the_worker->get_unnormalised_log_incremental_weights());
   }
   */

  // Comment for later...
  // IS: simulate
  // SMC w MCMC: t=1 simulate, t>1 loop resample-move, then weight, until stopping point reached with resample-move being the final step
  // PF/SMC: t=1 simulate, t>1 resample, sim prop, then weight, until stopping point reached with sim prop being final step
  // void ImportanceSampler::smc_simulate(SMCOutput* current_state)
  //{
  // this->output->add_proposed_particles(particles);

  // the_worker->simulate_and_weight();
  /*
   unsigned int number_of_points = algorithm["number_of_points"];

   List observed_data = model["observed_data"];

   // Do the initial importance sampling step.

   // Do the simulation.
   SEXP simulate_proposal_SEXP = algorithm["simulate_proposal"];
   SimulateImportanceSamplingProposalPtr simulate_proposal = load_simulate_distribution(simulate_proposal_SEXP);

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

  //}

  // void ImportanceSampler::smc_weight(SMCOutput* current_state)
  //{
  //   this->the_worker->weight();
  //   current_state->update_unnormalised_log_incremental_weights(this->get_unnormalised_log_incremental_weights());

  //  current_state->update_unnormalised_log_weights();
  //}
}
