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
    this->simulate_smc(simulation);
    this->evaluate_smc(simulation);
    simulation->normalise_and_resample_weights();
    return simulation;
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
    // set sequencer to have values from conditioned_on_parameters
    if (!this->sequencer_limit_is_fixed)
      this->sequencer.set_next_with_parameter(conditioned_on_parameters);
    else
      this->sequencer.reset();

    this->set_abc_scale(&current_state->back());

    // this->smcadaptive_given_smcfixed_weight(current_state, conditioned_on_parameters);
    // this->sequencer.find_next_target_bisection(current_state,this->index,conditioned_on_parameters);
    this->sequencer.find_next_target_bisection(current_state, this->index);
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

    current_state->log_likelihood = current_state->calculate_latest_log_normalising_constant_ratio();

    // if (this->sequencer_parameters!=NULL)
    //   current_state->back().schedule_parameters = *this->sequencer_parameters;
    current_state->back().schedule_parameters = this->sequencer.schedule_parameters;

    current_state->set_time_and_reset_start();

    current_state->terminate();
  }

  void ImportanceSampler::subsample_evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput *current_state,
                                                                                 const Parameters &conditioned_on_parameters)
  {
    // set sequencer to have values from conditioned_on_parameters
    if (!this->sequencer_limit_is_fixed)
      this->sequencer.set_next_with_parameter(conditioned_on_parameters);
    else
      this->sequencer.reset();

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

    current_state->log_likelihood = current_state->calculate_latest_log_normalising_constant_ratio();

    // if (this->sequencer_parameters!=NULL)
    //   current_state->back().schedule_parameters = *this->sequencer_parameters;
    current_state->back().schedule_parameters = this->sequencer.schedule_parameters;

    current_state->set_time_and_reset_start();

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
