#include <numeric>
#include <algorithm>
#include "smc_mcmc_move.h"
#include "single_point_move_output.h"
#include "adaptive_iterations_mcmc_termination.h"
#include "esjd_adaptor.h"
#include "filesystem.h"
#include <fstream>
#include "ensemble_factor_variables.h"
#include "chunk_parallel.h"
#include "progress_reporter.h"
#include <sstream>
#include <mutex>
#include <memory>
#include "factors.h"
#include "factor_variables.h"
#include "smc_output.h"
#include "cess_smc_criterion.h"
#include "exact_likelihood_estimator.h"
#include "annealed_likelihood_estimator.h"
#include "utils.h"
#include "parameter_particle_simulator.h"
#include "sequential_smc_worker.h"
#include "mcmc.h"
#include "move_output.h"
#include "standard_mcmc_output.h"
#include "custom_distribution_proposal_kernel.h"
#include "vector_factors.h"
#include "vector_index.h"

namespace ilike
{

// The factors in full_index that are not in used_index: those the SMC never
// evaluates because they cancel in its weights (typically the prior, when the
// prior is the initial proposal). Returns NULL if there are none. The caller
// owns the result.
static Index* make_cancelled_index(const Index* full_index,
                                   const Index* used_index)
{
  if ( (full_index==NULL) || (used_index==NULL) )
    return NULL;
  std::vector<size_t> cancelled;
  for (auto i=full_index->begin(); i!=full_index->end(); ++i)
  {
    if (std::find(used_index->begin(), used_index->end(), *i)==used_index->end())
      cancelled.push_back(*i);
  }
  if (cancelled.size()==0)
    return NULL;
  return new VectorIndex(cancelled);
}

SMCMCMCMove::SMCMCMCMove()
:SMC()
{
  this->mcmc.reset();
  this->index.reset();
  this->cancelled_index.reset();
}

// Multiple MCMC chains.
SMCMCMCMove::SMCMCMCMove(RandomNumberGenerator* rng_in,
                         size_t* seed_in,
                         Data* data_in,
                         const Parameters &algorithm_parameters,
                         size_t lag_in,
                         size_t lag_proposed_in,
                         MCMC* mcmc_in,
                         EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
                         EvaluateLogDistributionPtr evaluate_log_prior_in,
                         const std::vector<Parameters> &initial_points_in,
                         const arma::colvec &log_probabilities_of_initial_values_in,
                         bool transform_proposed_particles,
                         bool parallel_in,
                         size_t grain_size_in,
                         const std::string &results_name_in)
:SMC(rng_in,
     seed_in,
     data_in,
     algorithm_parameters,
     initial_points_in.size(),
     lag_in,
     lag_proposed_in,
     mcmc_in->get_proposals(),
     -1.0,
     false,
     true,
     true,
     transform_proposed_particles,
     results_name_in)
{
  this->chunk_size = grain_size_in;
  mcmc_in->set_proposal_parameters(&this->algorithm_parameters);
  
  std::vector<LikelihoodEstimator*> likelihood_estimators;
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
  
  this->particle_simulator.reset();
  
  this->initial_particles = initial_points_in;
  this->log_probabilities_of_initial_values = log_probabilities_of_initial_values_in;
  this->proposed_particles_inputted = true;
  
  if (initial_points_in.size()==0)
    Rcpp::stop("SMCMCMCMove::SMCMCMCMove - need more than one initialising point.");

  this->parallel = parallel_in;
  this->the_worker.reset(new SequentialSMCWorker(this));
  
  std::vector<double> schedule_in;
  schedule_in.push_back(0.0);
  schedule_in.push_back(1.0);
  std::string variable_in = "power";
  SMCCriterion* smc_criterion = new CESSSMCCriterion(-1.0);
  this->sequencer = Sequencer(this->the_worker.get(),
                              schedule_in,
                              variable_in,
                              25,
                              smc_criterion);
  
  this->mcmc.reset(mcmc_in);
  this->cancelled_index.reset();
  VectorIndex mcmc_index(indices); // set_index copies it
  this->mcmc->set_index(&mcmc_index);
  this->mcmc_at_last_step = true;
}

// Multiple MCMC chains.
SMCMCMCMove::SMCMCMCMove(RandomNumberGenerator* rng_in,
                         size_t* seed_in,
                         Data* data_in,
                         const Parameters &algorithm_parameters,
                         size_t number_of_particles_in,
                         size_t lag_in,
                         size_t lag_proposed_in,
                         MCMC* mcmc_in,
                         const std::vector<LikelihoodEstimator*> &likelihood_estimators_in,
                         IndependentProposalKernel* proposal_in,
                         Index* without_cancelled_index,
                         Index* full_index,
                         bool transform_proposed_particles,
                         bool parallel_in,
                         size_t grain_size_in,
                         const std::string &results_name_in)
:SMC(rng_in,
     seed_in,
     data_in,
     algorithm_parameters,
     number_of_particles_in,
     lag_in,
     lag_proposed_in,
     mcmc_in->get_proposals(),
     -1.0,
     true,
     true,
     true,
     transform_proposed_particles,
     results_name_in)
{
  this->chunk_size = grain_size_in;
  mcmc_in->set_proposal_parameters(&this->algorithm_parameters);
  proposal_in->set_proposal_parameters(&this->algorithm_parameters);
  
  this->index.reset(without_cancelled_index);
  
  this->factors.reset(new VectorFactors(likelihood_estimators_in));
  
  this->particle_simulator.reset(new ParameterParticleSimulator(proposal_in,
                                                            likelihood_estimators_in));
  
  this->proposed_particles_inputted = false;

  this->parallel = parallel_in;
  this->the_worker.reset(new SequentialSMCWorker(this));
  
  std::vector<double> schedule_in;
  schedule_in.push_back(0.0);
  schedule_in.push_back(1.0);
  std::string variable_in = "power";
  SMCCriterion* smc_criterion = new CESSSMCCriterion(-1.0);
  this->sequencer = Sequencer(this->the_worker.get(),
                              schedule_in,
                              variable_in,
                              25,
                              smc_criterion);
  
  this->mcmc.reset(mcmc_in);
  this->cancelled_index.reset();
  if (this->mcmc!=NULL)
    this->mcmc->set_index_if_null(full_index);
  this->cancelled_index.reset(make_cancelled_index(full_index, without_cancelled_index));
  this->mcmc_at_last_step = true;
}

SMCMCMCMove::SMCMCMCMove(RandomNumberGenerator* rng_in,
                         size_t* seed_in,
                         Data* data_in,
                         const Parameters &algorithm_parameters,
                         size_t lag_in,
                         size_t lag_proposed_in,
                         MCMC* mcmc_in,
                         const std::vector<LikelihoodEstimator*> &likelihood_estimators_in,
                         const std::vector<Parameters> &initial_points_in,
                         const arma::colvec &log_probabilities_of_initial_values_in,
                         Index* without_cancelled_index,
                         Index* full_index,
                         bool transform_proposed_particles,
                         bool parallel_in,
                         size_t grain_size_in,
                         const std::string &results_name_in)
:SMC(rng_in,
     seed_in,
     data_in,
     algorithm_parameters,
     initial_points_in.size(),
     lag_in,
     lag_proposed_in,
     mcmc_in->get_proposals(),
     -1.0,
     false,
     true,
     true,
     transform_proposed_particles,
     results_name_in)
{
  this->chunk_size = grain_size_in;
  mcmc_in->set_proposal_parameters(&this->algorithm_parameters);
  
  this->index.reset(without_cancelled_index);
  
  this->factors.reset(new VectorFactors(likelihood_estimators_in));
  
  this->particle_simulator.reset();
  
  this->initial_particles = initial_points_in;
  
  this->log_probabilities_of_initial_values = log_probabilities_of_initial_values_in;
  this->proposed_particles_inputted = true;
  
  if (initial_points_in.size()==0)
    Rcpp::stop("SMCMCMCMove::SMCMCMCMove - need more than one initialising point.");

  this->parallel = parallel_in;
  this->the_worker.reset(new SequentialSMCWorker(this));
  
  std::vector<double> schedule_in;
  schedule_in.push_back(0.0);
  schedule_in.push_back(1.0);
  std::string variable_in = "power";
  SMCCriterion* smc_criterion = new CESSSMCCriterion(-1.0);
  this->sequencer = Sequencer(this->the_worker.get(),
                              schedule_in,
                              variable_in,
                              100,
                              smc_criterion);
  
  this->mcmc.reset(mcmc_in);
  this->cancelled_index.reset();
  if (this->mcmc!=NULL)
    this->mcmc->set_index_if_null(full_index);
  this->cancelled_index.reset(make_cancelled_index(full_index, without_cancelled_index));
  this->mcmc_at_last_step = true;
}

SMCMCMCMove::SMCMCMCMove(RandomNumberGenerator* rng_in,
                         size_t* seed_in,
                         Data* data_in,
                         const Parameters &algorithm_parameters,
                         size_t number_of_particles_in,
                         size_t lag_in,
                         size_t lag_proposed_in,
                         MCMC* mcmc_in,
                         EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
                         EvaluateLogDistributionPtr evaluate_log_prior_in,
                         SimulateDistributionPtr simulate_proposal_in,
                         EvaluateLogDistributionPtr evaluate_log_proposal_in,
                         bool mcmc_at_last_step_in,
                         bool transform_proposed_particles,
                         bool parallel_in,
                         size_t grain_size_in,
                         const std::string &results_name_in)
:SMC(rng_in, seed_in, data_in, algorithm_parameters, number_of_particles_in, lag_in, lag_proposed_in, mcmc_in->get_proposals(), -1.0, true, true, true, transform_proposed_particles, results_name_in)
{
  this->chunk_size = grain_size_in;
  mcmc_in->set_proposal_parameters(&this->algorithm_parameters);
  
  IndependentProposalKernel* proposal = new CustomDistributionProposalKernel(simulate_proposal_in,
                                                                             evaluate_log_proposal_in);
  proposal->set_proposal_parameters(&this->algorithm_parameters);
  
  std::vector<LikelihoodEstimator*> likelihood_estimators;
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
  std::string variable_in = "power";
  SMCCriterion* smc_criterion = new CESSSMCCriterion(-1.0);
  this->sequencer = Sequencer(this->the_worker.get(),
                              schedule_in,
                              variable_in,
                              25,
                              smc_criterion);
  
  this->mcmc.reset(mcmc_in);
  this->cancelled_index.reset();
  VectorIndex mcmc_index(indices); // set_index copies it
  this->mcmc->set_index(&mcmc_index);
  this->mcmc_at_last_step = mcmc_at_last_step_in;
}

SMCMCMCMove::SMCMCMCMove(RandomNumberGenerator* rng_in,
                         size_t* seed_in,
                         Data* data_in,
                         const Parameters &algorithm_parameters,
                         size_t number_of_particles_in,
                         size_t lag_in,
                         size_t lag_proposed_in,
                         MCMC* mcmc_in,
                         double resampling_desired_ess_in,
                         double annealing_desired_cess_in,
                         size_t number_of_bisections_in,
                         EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
                         EvaluateLogDistributionPtr evaluate_log_prior_in,
                         SimulateDistributionPtr simulate_proposal_in,
                         EvaluateLogDistributionPtr evaluate_log_proposal_in,
                         bool mcmc_at_last_step_in,
                         bool transform_proposed_particles,
                         bool parallel_in,
                         size_t grain_size_in,
                         const std::string &results_name_in)
:SMC(rng_in, seed_in, data_in, algorithm_parameters, number_of_particles_in, lag_in, lag_proposed_in, mcmc_in->get_proposals(), resampling_desired_ess_in, true, true, true, transform_proposed_particles, results_name_in)
{
  this->chunk_size = grain_size_in;
  mcmc_in->set_proposal_parameters(&this->algorithm_parameters);
  
  std::string variable_in = "power";
  
  // Need to construct LikelihoodEstimator to read in to this constructor.
  IndependentProposalKernel* proposal = new CustomDistributionProposalKernel(simulate_proposal_in,
                                                                             evaluate_log_proposal_in);
  proposal->set_proposal_parameters(&this->algorithm_parameters);
  
  std::vector<LikelihoodEstimator*> likelihood_estimators;
  std::vector<size_t> indices;
  
  // Prior times likelihood.
  ExactLikelihoodEstimator* prior_times_likelihood = new ExactLikelihoodEstimator(rng_in,
                                                                                  seed_in,
                                                                                  data_in,
                                                                                  evaluate_log_prior_in,
                                                                                  evaluate_log_likelihood_in,
                                                                                  true);
  
  PowerFunctionPtr power = annealing_power;
  
  likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(rng_in,
                                                                  seed_in,
                                                                  data_in,
                                                                  prior_times_likelihood,
                                                                  power,
                                                                  variable_in,
                                                                  false));
  indices.push_back(0);
  
  // Proposal.
  ExactLikelihoodEstimator* proposal_for_evaluation = new ExactLikelihoodEstimator(rng_in,
                                                                                   seed_in,
                                                                                   data_in,
                                                                                   proposal,
                                                                                   true);
  
  PowerFunctionPtr second_power = annealing_one_minus_power;

  likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(rng_in,
                                                                  seed_in,
                                                                  data_in,
                                                                  proposal_for_evaluation,
                                                                  second_power,
                                                                  variable_in,
                                                                  false));
  
  indices.push_back(1);
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
  SMCCriterion* smc_criterion = new CESSSMCCriterion(annealing_desired_cess_in);
  this->sequencer = Sequencer(this->the_worker.get(),
                              schedule_in,
                              variable_in,
                              number_of_bisections_in,
                              smc_criterion);
  
  this->mcmc.reset(mcmc_in);
  this->cancelled_index.reset();
  VectorIndex mcmc_index(indices); // set_index copies it
  this->mcmc->set_index(&mcmc_index);
  this->mcmc_at_last_step = mcmc_at_last_step_in;
}

SMCMCMCMove::SMCMCMCMove(RandomNumberGenerator* rng_in,
                         size_t* seed_in,
                         Data* data_in,
                         const Parameters &algorithm_parameters,
                         size_t number_of_particles_in,
                         size_t lag_in,
                         size_t lag_proposed_in,
                         MCMC* mcmc_in,
                         double resampling_desired_ess_in,
                         double annealing_desired_cess_in,
                         size_t number_of_bisections_in,
                         const std::vector<double> &temperatures_in,
                         EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
                         EvaluateLogDistributionPtr evaluate_log_prior_in,
                         SimulateDistributionPtr simulate_proposal_in,
                         EvaluateLogDistributionPtr evaluate_log_proposal_in,
                         bool mcmc_at_last_step_in,
                         bool transform_proposed_particles,
                         bool parallel_in,
                         size_t grain_size_in,
                         const std::string &results_name_in)
:SMC(rng_in, seed_in, data_in, algorithm_parameters, number_of_particles_in, lag_in, lag_proposed_in, mcmc_in->get_proposals(), resampling_desired_ess_in, true, true, true, transform_proposed_particles, results_name_in)
{
  this->chunk_size = grain_size_in;
  mcmc_in->set_proposal_parameters(&this->algorithm_parameters);
  std::string variable_in = "power";
  // Need to construct LikelihoodEstimator to read in to this constructor.
  IndependentProposalKernel* proposal = new CustomDistributionProposalKernel(simulate_proposal_in,
                                                                             evaluate_log_proposal_in);
  proposal->set_proposal_parameters(&this->algorithm_parameters);
  
  std::vector<LikelihoodEstimator*> likelihood_estimators;
  std::vector<size_t> indices;
  
  // Prior times likelihood.
  ExactLikelihoodEstimator* prior_times_likelihood = new ExactLikelihoodEstimator(rng_in,
                                                                                  seed_in,
                                                                                  data_in,
                                                                                  evaluate_log_prior_in,
                                                                                  evaluate_log_likelihood_in,
                                                                                  true);
  
  PowerFunctionPtr power = annealing_power;
  
  likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(rng_in,
                                                                  seed_in,
                                                                  data_in,
                                                                  prior_times_likelihood,
                                                                  power,
                                                                  variable_in,
                                                                  false));
  indices.push_back(0);
  
  // Proposal.
  ExactLikelihoodEstimator* proposal_for_evaluation = new ExactLikelihoodEstimator(rng_in,
                                                                                   seed_in,
                                                                                   data_in,
                                                                                   proposal,
                                                                                   true);
  
  PowerFunctionPtr second_power = annealing_one_minus_power;

  likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(rng_in,
                                                                  seed_in,
                                                                  data_in,
                                                                  proposal_for_evaluation,
                                                                  second_power,
                                                                  variable_in,
                                                                  false));
  
  indices.push_back(1);
  this->index.reset(new VectorIndex(indices));
  
  this->factors.reset(new VectorFactors(likelihood_estimators));

  // Need to construct LikelihoodEstimator to read in to this constructor.
  this->particle_simulator.reset(new ParameterParticleSimulator(proposal,
                                                            likelihood_estimators));
  
  this->parallel = parallel_in;
  this->the_worker.reset(new SequentialSMCWorker(this));
  
  std::vector<double> schedule_in;
  schedule_in.reserve(temperatures_in.size()+1);
  schedule_in.push_back(0.0);
  for (size_t i=0; i<temperatures_in.size(); ++i)
  {
    schedule_in.push_back(temperatures_in[i]);
  }
  
  SMCCriterion* smc_criterion = new CESSSMCCriterion(annealing_desired_cess_in);
  this->sequencer = Sequencer(this->the_worker.get(),
                              schedule_in,
                              variable_in,
                              number_of_bisections_in,
                              smc_criterion);
  
  this->mcmc.reset(mcmc_in);
  this->cancelled_index.reset();
  VectorIndex mcmc_index(indices); // set_index copies it
  this->mcmc->set_index(&mcmc_index);
  this->mcmc_at_last_step = mcmc_at_last_step_in;
}

SMCMCMCMove::SMCMCMCMove(RandomNumberGenerator* rng_in,
                         size_t* seed_in,
                         Data* data_in,
                         const Parameters &algorithm_parameters,
                         size_t number_of_particles_in,
                         size_t lag_in,
                         size_t lag_proposed_in,
                         MCMC* mcmc_in,
                         double resampling_desired_ess_in,
                         double annealing_desired_cess_in,
                         size_t number_of_bisections_in,
                         SMCTermination* termination_in,
                         const std::string &sequence_variable_in,
                         const std::vector<double> &schedule_in,
                         const std::vector<LikelihoodEstimator*> &likelihood_estimators_in,
                         IndependentProposalKernel* proposal_in,
                         bool proposal_is_evaluated_in,
                         bool smcfixed_flag_in,
                         bool sequencer_limit_is_fixed_in,
                         bool mcmc_at_last_step_in,
                         bool transform_proposed_particles,
                         bool parallel_in,
                         size_t grain_size_in,
                         const std::string &results_name_in)
:SMC(rng_in, seed_in, data_in, algorithm_parameters, number_of_particles_in, lag_in, lag_proposed_in, mcmc_in->get_proposals(), resampling_desired_ess_in, proposal_is_evaluated_in, smcfixed_flag_in, sequencer_limit_is_fixed_in, transform_proposed_particles, results_name_in)
{
  this->chunk_size = grain_size_in;
  
  mcmc_in->set_proposal_parameters(&this->algorithm_parameters);
  proposal_in->set_proposal_parameters(&this->algorithm_parameters);
  
  std::vector<size_t> indices;
  indices.reserve(likelihood_estimators_in.size());
  for (size_t i=0; i<likelihood_estimators_in.size(); ++i)
    indices.push_back(i);
  this->index.reset(new VectorIndex(indices));
  
  this->factors.reset(new VectorFactors(likelihood_estimators_in));

  // Need to construct LikelihoodEstimator to read in to this constructor.
  this->particle_simulator.reset(new ParameterParticleSimulator(proposal_in,
                                                            likelihood_estimators_in));
  
  this->parallel = parallel_in;
  this->the_worker.reset(new SequentialSMCWorker(this));
  
  SMCCriterion* smc_criterion = new CESSSMCCriterion(annealing_desired_cess_in);
  this->sequencer = Sequencer(this->the_worker.get(),
                              schedule_in,
                              sequence_variable_in,
                              number_of_bisections_in,
                              smc_criterion,
                              termination_in);
  
  this->mcmc.reset(mcmc_in);
  this->cancelled_index.reset();
  if (this->mcmc!=NULL)
  {
    VectorIndex mcmc_index(indices); // set_index_if_null copies it
    this->mcmc->set_index_if_null(&mcmc_index);
  }
  this->mcmc_at_last_step = mcmc_at_last_step_in;
}

SMCMCMCMove::SMCMCMCMove(RandomNumberGenerator* rng_in,
                         size_t* seed_in,
                         Data* data_in,
                         const Parameters &algorithm_parameters,
                         size_t number_of_particles_in,
                         size_t lag_in,
                         size_t lag_proposed_in,
                         MCMC* mcmc_in,
                         SMCCriterion* adaptive_resampling_in,
                         SMCCriterion* adaptive_target_in,
                         size_t number_of_bisections_in,
                         SMCTermination* termination_in,
                         const std::vector<std::string> &sequence_variables_in,
                         const std::vector<std::vector<double>> &schedules_in,
                         const std::vector<LikelihoodEstimator*> &likelihood_estimators_in,
                         IndependentProposalKernel* proposal_in,
                         Index* without_cancelled_index,
                         Index* full_index,
                         bool proposal_is_evaluated_in,
                         bool smcfixed_flag_in,
                         bool sequencer_limit_is_fixed_in,
                         bool mcmc_at_last_step_in,
                         bool transform_proposed_particles,
                         bool parallel_in,
                         size_t grain_size_in,
                         const std::string &results_name_in)
:SMC(rng_in, seed_in, data_in, algorithm_parameters, number_of_particles_in, lag_in, lag_proposed_in, mcmc_in->get_proposals(), adaptive_resampling_in, proposal_is_evaluated_in, smcfixed_flag_in, sequencer_limit_is_fixed_in, transform_proposed_particles, results_name_in)
{
  this->chunk_size = grain_size_in;
  
  mcmc_in->set_proposal_parameters(&this->algorithm_parameters);
  proposal_in->set_proposal_parameters(&this->algorithm_parameters);
  
  this->index.reset(without_cancelled_index);
  
  this->factors.reset(new VectorFactors(likelihood_estimators_in));

  // Need to construct LikelihoodEstimator to read in to this constructor.
  this->particle_simulator.reset(new ParameterParticleSimulator(proposal_in,
                                                            likelihood_estimators_in));

  this->parallel = parallel_in;
  this->the_worker.reset(new SequentialSMCWorker(this));
  
  SMCCriterion* smc_criterion = adaptive_target_in;
  this->sequencer = Sequencer(this->the_worker.get(),
                              schedules_in,
                              sequence_variables_in,
                              number_of_bisections_in,
                              smc_criterion,
                              termination_in);
  
  this->mcmc.reset(mcmc_in);
  this->cancelled_index.reset();
  if (this->mcmc!=NULL)
    this->mcmc->set_index_if_null(full_index);
  this->cancelled_index.reset(make_cancelled_index(full_index, without_cancelled_index));
  
  this->mcmc_at_last_step = mcmc_at_last_step_in;
}

//Copy constructor for the SMCMCMCMove class.
SMCMCMCMove::SMCMCMCMove(const SMCMCMCMove &another)
:SMC(another)
{
  this->make_copy(another);
}

//Destructor for the SMCMCMCMove class.
SMCMCMCMove::~SMCMCMCMove()
{
  this->mcmc.reset();
  
  this->index.reset();
  
  this->cancelled_index.reset();
}

void SMCMCMCMove::operator=(const SMCMCMCMove &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  this->mcmc.reset();
  
  this->index.reset();
  
  this->cancelled_index.reset();
  
  SMC::operator=(another);
  this->make_copy(another);
}

SMC* SMCMCMCMove::smc_duplicate(void) const
{
  return( new SMCMCMCMove(*this));
}

LikelihoodEstimator* SMCMCMCMove::duplicate() const
{
  return( new SMCMCMCMove(*this));
}

void SMCMCMCMove::make_copy(const SMCMCMCMove &another)
{
  if (another.mcmc!=NULL)
    this->mcmc.reset(another.mcmc->mcmc_duplicate());
  else
    this->mcmc.reset();
  
  if (another.index!=NULL)
    this->index.reset(another.index->duplicate());
  else
    this->index.reset();
  
  if (another.cancelled_index!=NULL)
    this->cancelled_index.reset(another.cancelled_index->duplicate());
  else
    this->cancelled_index.reset();
  
  this->mcmc_at_last_step = another.mcmc_at_last_step;
  this->move_only_after_resampling = another.move_only_after_resampling;
  this->follows_caller_target = another.follows_caller_target;
  this->followed_variable = another.followed_variable;
  this->followed_targets = another.followed_targets;
  // (the tuned kernels are this copy's own: found again)
  this->tuners.clear();
  this->tuners_found = false;
  this->particle_tunings.clear();
}

// esjd (C211): the tuned kernels' values for this move step, and each
// particle's rows; with a summary of the values in results/esjd.csv
void SMCMCMCMove::prepare_tuning(SMCOutput* current_state)
{
  if (!this->tuners_found)
  {
    this->tuners_found = true;
    this->tuners.clear();
    for (const ProposalKernel* kernel : this->mcmc->get_proposals())
    {
      SMCAdaptor* adaptor = kernel->get_smc_adaptor();
      std::vector<SMCAdaptor*> adaptors;
      if (CompositeSMCAdaptor* composite = dynamic_cast<CompositeSMCAdaptor*>(adaptor))
        adaptors = composite->get_parts();
      else if (adaptor != NULL)
        adaptors.push_back(adaptor);
      for (SMCAdaptor* a : adaptors)
        if (ESJDAdaptor* tuner = dynamic_cast<ESJDAdaptor*>(a))
          this->tuners.push_back(std::make_pair(kernel->get_instance_index(), tuner));
    }
  }
  if (this->tuners.empty())
    return;
  if (this->chains_only)
    for (auto &t : this->tuners)
      if (t.second->is_choice())
        Rcpp::stop("choose([...]) is used in SMC_with_MCMC, where the particles show which proposal moves them best "
                   "(in MCMC each chain would keep the proposal drawn at its start).");
  if (this->waste_free_chains > 0)
    Rcpp::stop("smc_adapt = esjd(…) is not available with waste-free SMC.");
  if (!this->can_move_flat() && (this->rich_chunk_mode() == 0))
    Rcpp::stop("smc_adapt = esjd(…) is not available for this model's SMC moves (the factors must be exact, "
               "annealed, or support parallel chunks or copies).");

  size_t N = this->number_of_particles;
  // its own stream from the run's seed (the moves' streams are unchanged)
  RandomNumberGenerator tuning_rng(*this->the_worker->get_rng());
  tuning_rng.seed(this->the_worker->get_seed() ^ 0x5DEECE66DULL);
  for (auto &t : this->tuners)
    t.second->update(tuning_rng, current_state, N);
  // pretune (C223): trial moves from each particle's starting state (as the
  // move starts it: its ancestor, restored on a flat generation, with the
  // cancelled factors evaluated)
  bool trials = false;
  for (auto &t : this->tuners)
    trials = trials || t.second->needs_trial_moves();
  if (trials)
  {
    Particles &generation = current_state->back();
    bool from_flat = generation.is_flat();
    const Parameters &restore_parameters =
      (from_flat && generation.adaptive_summaries_used) ? generation.last_trial_schedule
                                                         : this->sequencer.schedule_parameters;
    auto start = [&](size_t i) -> Particle {
      size_t ancestor = generation.ancestor_variables[i];
      Particle particle;
      if (from_flat)
        this->restore_particle(particle, generation, ancestor, restore_parameters);
      else
        particle = generation[ancestor]->back();
      if (this->cancelled_index != NULL)
        particle.evaluate_likelihoods(this->cancelled_index.get());
      return particle;
    };
    for (auto &t : this->tuners)
      if (t.second->needs_trial_moves())
        t.second->trial_moves(tuning_rng, N, start);
  }
  if (this->particle_tunings.size() != N)
    this->particle_tunings.resize(N);
  for (size_t i = 0; i < N; ++i)
  {
    ParticleTuning &tuning = this->particle_tunings[i];
    tuning.rows.clear();
    for (auto &t : this->tuners)
    {
      ParticleTuning::Row row;
      row.instance = t.first;
      row.values = t.second->column(i);
      row.adaptor = t.second;
      row.choice = t.second->is_choice();
      tuning.rows.push_back(row);
    }
  }

  // the values' summary, one row per move step and tuned kernel
  if (!this->results_name.empty())
  {
    // (beside the SMC's output.h5, in <results>/ilike_smc)
    std::string directory = this->results_name + "_smc";
    if (!directory_exists(directory))
      make_directory(directory);
    std::string path = directory + "/esjd.csv";
    bool exists = std::ifstream(path).good();
    std::ofstream file(path, std::ios::app);
    if (file.good())
    {
      // (long format: one row per statistic, as tuners differ in their
      // statistics; load_smc_tuning widens it, C224)
      if (!exists)
        file << "step,tuner,variable,statistic,value\n";
      for (size_t k = 0; k < this->tuners.size(); ++k)
      {
        std::vector<std::string> names;
        std::vector<double> numbers;
        this->tuners[k].second->summary(names, numbers);
        for (size_t j = 0; j < names.size(); ++j)
          file << this->tuners[k].second->get_updates() << "," << (k + 1) << ","
               << stacked_display_name(this->tuners[k].second->get_variable()) << "," << names[j] << "," << numbers[j] << "\n";
      }
    }
  }
}

SMCOutput* SMCMCMCMove::specific_run()
{
  SMCOutput* simulation = this->initialise_smc();
  // progress (dev/plan.md 3.11): MCMC() counts iterations over the chains;
  // SMC shows the position along its sequence of targets (percent)
  double progress_total = 100.0;
  std::string progress_name = "SMC";
  if (this->progress_counts_iterations)
  {
    size_t iterations = (this->mcmc != NULL) ? this->mcmc->known_number_of_iterations() : 0;
    progress_total = (iterations > 0) ? double(iterations)*double(this->number_of_particles) : -1.0;
    progress_name = "MCMC";
  }
  ProgressScope progress(this->show_progress, progress_total, progress_name);
  this->simulate_smc(simulation);
  this->evaluate_smc(simulation);
  simulation->normalise_and_resample_weights();
  return simulation;
}

SMCOutput* SMCMCMCMove::specific_initialise_smc()
{
  SMCOutput* output = new SMCOutput(this, this->lag, this->lag_proposed, this->results_name);
  return output;
}

int SMCMCMCMove::rich_chunk_mode() const
{
  if ((dynamic_cast<const SequentialSMCWorker*>(this->the_worker.get()) == NULL) ||
      (this->mcmc == NULL) || (this->factors == NULL))
    return 0;
  if (this->factors->supports_parallel_chunks())
    return 1;
  if (this->factors->supports_parallel_copies())
    return 2;
  return 0;
}

void SMCMCMCMove::move_rich_in_chunks(Particles* next_particles,
                                      Particles* current_particles,
                                      bool copies)
{
  // SequentialSMCWorker::specific_move and the seed bookkeeping of
  // SMCWorker::move, grain_size particles at a time, each chunk with its own
  // stream, each particle's MoveOutput put in its own slot. A kernel that
  // adapts is copied for every chain (C98), which touches a static counter
  // of the kernels: such runs are chunked but not parallel.
  SMCWorker* worker = this->the_worker.get();
  size_t N = this->number_of_particles;
  size_t B = ((this->chunk_size == 0) || (this->chunk_size > N)) ? N : this->chunk_size;
  size_t K = (N + B - 1)/B;
  bool in_parallel = this->parallel && !this->mcmc->adapts();
  if (this->parallel && !in_parallel)
    this->parallel_refused = true;
  next_particles->particles.clear();
  next_particles->particles.resize(N);

  RandomNumberGenerator first_stream(*worker->get_rng());
  first_stream.seed(worker->get_seed());
  std::vector<RandomNumberGenerator> streams = chunk_streams(first_stream, K);

  run_chunk_groups(0, K, in_parallel, this->number_of_threads,
                   [&](size_t group, size_t first_chunk, size_t end_chunk)
                   {
                     // the group's copy of the factors, drawing from its own
                     // generator and seed counter
                     std::unique_ptr<Factors> own_factors;
                     RandomNumberGenerator group_rng;
                     size_t group_seed = 0;
                     const Factors* factors_for_group = this->factors.get();
                     if (copies)
                     {
                       own_factors.reset(this->factors->duplicate());
                       own_factors->redirect_random_numbers(&group_rng, &group_seed);
                       factors_for_group = own_factors.get();
                     }
                     for (size_t k = first_chunk; k < end_chunk; ++k)
                     {
                       size_t first = k*B;
                       size_t n = std::min(B, N - first);
                       RandomNumberGenerator local_rng(streams[k]);
                       for (size_t i = first; i < first + n; ++i)
                       {
                         const Particle &ancestor = (*current_particles)[current_particles->ancestor_variables[i]]->back();
                         MoveOutput* moved;
                         bool tuned = !this->tuners.empty();
                         if (copies)
                         {
                           group_seed = size_t(local_rng());
                           Particle start = ancestor;
                           if (tuned)
                             start.tuning = &this->particle_tunings[i];
                           if ((start.factor_variables != NULL) && !start.factor_variables->repoint_to(factors_for_group))
                             Rcpp::stop("SMCMCMCMove::move_rich_in_chunks: the factors could not be repointed.");
                           moved = this->move(local_rng, start);
                           // what is kept outlives the group's copy
                           if (!moved->repoint_to(this->factors.get()))
                             Rcpp::stop("SMCMCMCMove::move_rich_in_chunks: the factors could not be repointed.");
                         }
                         else if (tuned)
                         {
                           Particle start = ancestor;
                           start.tuning = &this->particle_tunings[i];
                           moved = this->move(local_rng, start);
                         }
                         else
                         {
                           moved = this->move(local_rng, ancestor);
                         }
                         if (tuned)
                         {
                           for (auto &t : this->tuners)
                             t.second->score(i, ancestor.parameters, moved->back().parameters, this->particle_tunings[i]);
                           moved->back().tuning = NULL;
                         }
                         moved->back().previous_self = &ancestor;
                         next_particles->particles[i].reset(moved);
                       }
                     }
                   },
                   [&](size_t chunks_done)
                   {
                     if (this->progress_counts_iterations)
                     {
                       size_t iterations = this->mcmc->known_number_of_iterations();
                       ProgressReporter::instance().update(double(std::min(N, chunks_done*B))*double(iterations));
                     }
                   });

  // as SMCWorker::move
  worker->set_seed(worker->get_seed() + N);
  this->rng->seed(worker->get_seed());
  next_particles->simulate_resampling_variables(*this->rng);
  worker->set_seed(worker->get_seed() + N);
}

bool SMCMCMCMove::can_simulate_flat() const
{
  // the first generation can be simulated flat when the moves can start
  // from it (can_move_flat), its smcfixed parts use no random numbers and
  // give summaries, and the particles are simulated (not given)
  return !this->proposed_particles_inputted && (this->particle_simulator != NULL) &&
    (this->number_of_particles > 0) && this->can_move_flat() &&
    this->factors->supports_flat_generations(this->index.get());
}

bool SMCMCMCMove::can_move_flat() const
{
  // Conditions (dev/plan.md 2.6a, G3): the moves read the cached factor
  // values (so the cancelled factors are evaluated before each run, as in
  // move()); one kernel shared by all chains (no adaptation); a run length
  // fixed in advance (iterations(n): no termination reading the chain); and
  // every factor restorable from its summary (the SMC's factors) or
  // reinitialisable (the cancelled ones).
  if (this->force_rich_generations || this->move_only_after_resampling)
    return false;
  if (dynamic_cast<const SequentialSMCWorker*>(this->the_worker.get()) == NULL)
    return false;
  if ((this->mcmc == NULL) || !this->mcmc->reads_cached_factors() || this->mcmc->adapts())
    return false;
  if (this->mcmc->known_number_of_iterations() == 0)
    return false;
  // moves with per-iteration diagnostics (nuts, C212) record them in each
  // chain's output, so they run on the rich path (a nuts iteration costs many
  // gradients: the flat path's saving is negligible there)
  std::vector<std::string> diagnostics;
  this->mcmc->diagnostic_names(diagnostics);
  if (!diagnostics.empty())
    return false;
  if ((this->factors == NULL) || !this->factors->supports_flat_moves(this->index.get(), this->cancelled_index.get()))
    return false;
  return true;
}

void SMCMCMCMove::restore_particle(Particle &particle,
                                   const Particles &generation,
                                   size_t row,
                                   const Parameters &restore_parameters) const
{
  // What a copy of the rich ancestor would be, for everything the moves
  // read: its variables (the SMC sequence variables stay shared with the
  // sequencer, as in every particle of this SMC), its factor outputs as its
  // last evaluation left them, and a new proposal store (gradients and
  // transforms are recomputed from the variables when needed).
  if (particle.parameters.vector_size() == 0)
    particle.parameters.merge_with_fixed(this->sequencer.schedule_parameters);
  // a generation made by moves keeps its chains' final states; one simulated
  // from the proposal keeps its points (one per particle)
  const SampleStore &points = (generation.flat_final_points.size() > 0) ? generation.flat_final_points
                                                                        : generation.flat_points;
  points.assign_into(row, particle.parameters);
  if (particle.factor_variables == NULL)
    particle.factor_variables.reset(this->factors->simulate_factor_variables(particle.parameters));
  particle.factor_variables->set_particle(&particle);
  if (!particle.factor_variables->restore_from_summary(this->index.get(),
                                                      this->cancelled_index.get(),
                                                      generation.adaptive_summaries.colptr(row),
                                                      restore_parameters))
    Rcpp::stop("SMCMCMCMove::restore_particle: the factors of a particle could not be restored.");
  if (particle.ensemble_factor_variables != NULL)
  {
    particle.ensemble_factor_variables.reset();
    particle.ensemble_factor_variables.reset();
  }
  particle.current_proposal_store.clear();
  particle.simulate_proposal_variables(&this->proposals_to_transform_for,
                                       &this->proposals_to_find_gradient_for);
  particle.accepted_outputs.clear();
  particle.previous_self = NULL;
  particle.previous_target_evaluated = generation.adaptive_previous_targets[row];
  particle.target_evaluated = generation.adaptive_previous_targets[row];
  particle.subsample_previous_target_evaluated = 0.0;
  particle.subsample_target_evaluated = 0.0;
  particle.previous_ensemble_target_evaluated = 0.0;
  particle.subsample_previous_ensemble_target_evaluated = 0.0;
  particle.ensemble_target_evaluated = 0.0;
  particle.subsample_ensemble_target_evaluated = 0.0;
}

void SMCMCMCMove::move_waste_free(Particles* next_particles,
                                  Particles* current_particles)
{
  // Waste-free SMC (Dau & Chopin 2020, Algorithm 2; dev/plan.md 2.7). The N
  // weighted points of the current generation give M ancestors, by
  // stratified resampling with the generation's first M resampling
  // variables; chain m starts at ancestor a_m and makes P - 1 = N/M - 1 MCMC
  // moves (the MCMC termination is not used); all its P states, the
  // ancestor included, are particles of the next generation, each with the
  // values a moved particle carries (its variables, factor records, adaptive
  // summary and its target over the SMC's index, the next weights' previous
  // target). The ancestor of particle m*P + s is a_m. Chunks are groups of
  // chains, each with its own stream; in parallel when asked.
  SMCWorker* worker = this->the_worker.get();
  size_t N = this->number_of_particles;
  size_t M = this->waste_free_chains;
  if ((M == 0) || (M > N) || ((N % M) != 0))
    Rcpp::stop("SMC_with_MCMC: waste_free_chains must divide the number of particles.");
  size_t P = N/M;
  size_t B = ((this->chunk_size == 0) || (this->chunk_size > M)) ? M : this->chunk_size;
  size_t K = (M + B - 1)/B;
  bool from_flat = current_particles->is_flat();
  const Parameters &restore_parameters =
    (from_flat && current_particles->adaptive_summaries_used) ? current_particles->last_trial_schedule
                                                               : this->sequencer.schedule_parameters;

  // M ancestors from the N weighted points (stratified), and the ancestor of
  // every particle of the next generation (for the output)
  arma::colvec weights = exp(current_particles->normalised_log_weights - log_sum_exp(current_particles->normalised_log_weights));
  arma::colvec cumulative = cumsum(weights/arma::accu(weights));
  std::vector<size_t> chain_ancestors(M);
  size_t j = 0;
  for (size_t m = 0; m < M; ++m)
  {
    double v = (current_particles->resampling_variables[m] + double(m))/double(M);
    while ((j < N - 1) && (cumulative[j] < v))
      j = j + 1;
    chain_ancestors[m] = j;
  }
  current_particles->ancestor_variables.assign(N, 0);
  for (size_t m = 0; m < M; ++m)
    for (size_t s = 0; s < P; ++s)
      current_particles->ancestor_variables[m*P + s] = chain_ancestors[m];
  current_particles->resampled_flag = true;
  // every step resamples: the new particles start with equal weights
  // (add_particles copied the current weights if the resampling criterion
  // had not resampled)
  next_particles->previous_normalised_log_weights.set_size(N);
  next_particles->previous_normalised_log_weights.fill(-log(double(N)));

  next_particles->flat = true;
  next_particles->flat_size = N;
  next_particles->adaptive_previous_targets.set_size(N);

  bool keep_records = this->write_factor_output;
  size_t number_of_records = 0;
  int summary_size = -1;
  std::once_flag setup;

  RandomNumberGenerator first_stream(*worker->get_rng());
  first_stream.seed(worker->get_seed());
  std::vector<RandomNumberGenerator> streams = chunk_streams(first_stream, K);

  auto keep_state = [&](Particle &state, size_t row, std::vector<double> &records)
  {
    const FactorVariables* factor_variables = state.factor_variables.get();
    int size = (factor_variables == NULL) ? -1 : factor_variables->adaptive_summary_size(this->index.get());
    std::call_once(setup, [&]()
    {
      if (keep_records && (state.factor_variables != NULL))
      {
        if (!state.factor_variables->has_compact_values())
          Rcpp::stop("SMCMCMCMove::move_waste_free: a factor has no compact record.");
        number_of_records = state.factor_variables->number_of_compact_values();
        next_particles->flat_points.set_number_of_extra_columns(number_of_records);
        next_particles->flat_record_writer.reset(state.factor_variables->duplicate());
        next_particles->flat_record_writer->set_particle(NULL);
      }
      next_particles->flat_points.prepare(state.parameters, N);
      if (size < 0)
        Rcpp::stop("SMCMCMCMove::move_waste_free: a factor has no adaptive summary.");
      summary_size = size;
      next_particles->adaptive_summaries.set_size(size_t(summary_size), N);
      next_particles->flat_summary_evaluator.reset(factor_variables->duplicate());
      next_particles->flat_summary_evaluator->set_particle(NULL);
      next_particles->adaptive_summary_evaluator = next_particles->flat_summary_evaluator.get();
      next_particles->adaptive_summary_index = this->index.get();
      next_particles->adaptive_summary_state = 1;
    });
    if (size != summary_size)
      Rcpp::stop("SMCMCMCMove::move_waste_free: the particles' adaptive summaries differ in size.");
    if (number_of_records > 0)
    {
      records.resize(number_of_records);
      state.factor_variables->compact_values(records.data());
      next_particles->flat_points.set_point(row, state.parameters, records.data());
    }
    else
    {
      next_particles->flat_points.set_point(row, state.parameters);
    }
    if (summary_size > 0)
      factor_variables->adaptive_summary(this->index.get(), next_particles->adaptive_summaries.colptr(row));
    next_particles->adaptive_previous_targets[row] = state.get_cached_likelihoods(this->index.get());
  };

  run_chunk_groups(0, K, this->parallel, this->number_of_threads,
                   [&](size_t group, size_t first_chunk, size_t end_chunk)
                   {
                     Particle current;
                     Particle scratch;
                     std::vector<double> records;
                     for (size_t k = first_chunk; k < end_chunk; ++k)
                     {
                       size_t first = k*B;
                       size_t n = std::min(B, M - first);
                       RandomNumberGenerator local_rng(streams[k]);
                       for (size_t m = first; m < first + n; ++m)
                       {
                         size_t ancestor = chain_ancestors[m];
                         if (from_flat)
                           this->restore_particle(current, *current_particles, ancestor, restore_parameters);
                         else
                           current = (*current_particles)[ancestor]->back();
                         if (this->cancelled_index != NULL)
                           current.evaluate_likelihoods(this->cancelled_index.get());
                         keep_state(current, m*P, records);
                         for (size_t s = 1; s < P; ++s)
                         {
                           this->mcmc->move_in_place(local_rng, current, scratch);
                           current.tell_factors_to_forget_they_were_already_written_to_file();
                           keep_state(current, m*P + s, records);
                         }
                       }
                     }
                   });

  // as SMCWorker::move
  worker->set_seed(worker->get_seed() + N);
  this->rng->seed(worker->get_seed());
  next_particles->simulate_resampling_variables(*this->rng);
  worker->set_seed(worker->get_seed() + N);
}

void SMCMCMCMove::move_flat(Particles* next_particles,
                            Particles* current_particles)
{
  // The steps of SequentialSMCWorker::specific_move with this->move() and
  // MCMC::run for each particle, then the seed bookkeeping of SMCWorker::move.
  // The chunks run in order, or in parallel (G4, G5), each with its own
  // stream and writing only its own particles' rows: a chain of
  // `iterations` moves stores iterations/thin states, so particle i's states
  // are rows [i*m, (i+1)*m) of the chain store.
  SMCWorker* worker = this->the_worker.get();
  size_t N = this->number_of_particles;
  size_t B = ((this->chunk_size == 0) || (this->chunk_size > N)) ? N : this->chunk_size;
  size_t K = (N + B - 1)/B;
  size_t iterations = this->mcmc->known_number_of_iterations();
  size_t thin = std::max<size_t>(this->mcmc->get_thin(), 1);
  size_t m = iterations/thin;
  bool from_flat = current_particles->is_flat();
  // the factor values a rich generation's particles held after the target
  // was chosen: those of the last trial target (SMC::finish_weighting_from_summaries)
  const Parameters &restore_parameters =
    (from_flat && current_particles->adaptive_summaries_used) ? current_particles->last_trial_schedule
                                                               : this->sequencer.schedule_parameters;

  next_particles->flat = true;
  next_particles->flat_size = N;
  next_particles->flat_output_lengths.assign(N, double(m));
  next_particles->adaptive_previous_targets.set_size(N);

  bool keep_records = this->write_factor_output;
  size_t number_of_records = 0;
  int summary_size = -1;
  // the first state to be stored, and the first particle to finish, fix the
  // layouts (every particle gives the same; set_point checks it)
  std::once_flag chain_setup;
  std::once_flag final_setup;

  RandomNumberGenerator first_stream(*worker->get_rng());
  first_stream.seed(worker->get_seed());
  std::vector<RandomNumberGenerator> streams = chunk_streams(first_stream, K);

  auto run_chunk = [&](size_t k, Particle &current, Particle &scratch, std::vector<double> &records)
  {
    size_t first = k*B;
    size_t n = std::min(B, N - first);
    RandomNumberGenerator local_rng(streams[k]);

    for (size_t i = first; i < first + n; ++i)
    {
      size_t ancestor = current_particles->ancestor_variables[i];
      if (from_flat)
        this->restore_particle(current, *current_particles, ancestor, restore_parameters);
      else
        current = (*current_particles)[ancestor]->back();

      // as move(): the cancelled factors, then the run
      if (this->cancelled_index != NULL)
        current.evaluate_likelihoods(this->cancelled_index.get());
      // (esjd, C211: the particle's own tuning, and where it starts)
      bool tuned = !this->tuners.empty();
      Parameters before;
      if (tuned)
      {
        current.tuning = &this->particle_tunings[i];
        before = current.parameters;
      }

      size_t stored = 0;
      for (size_t iteration = 1; iteration <= iterations; ++iteration)
      {
        this->mcmc->move_in_place(local_rng, current, scratch);
        current.tell_factors_to_forget_they_were_already_written_to_file();
        if ((iteration % 1024) == 0)
        {
          check_interrupt(); // in long chains (dev/plan.md 3.11)
          if (this->progress_counts_iterations)
            ProgressReporter::instance().update(double(i)*double(iterations) + double(iteration));
        }
        if ((thin > 1) && ((iteration % thin) != 0))
          continue;
        // as StandardMCMCOutput::record
        std::call_once(chain_setup, [&]()
        {
          if (keep_records && (current.factor_variables != NULL))
          {
            if (!current.factor_variables->has_compact_values())
              Rcpp::stop("SMCMCMCMove::move_flat: a factor has no compact record.");
            number_of_records = current.factor_variables->number_of_compact_values();
            next_particles->flat_points.set_number_of_extra_columns(number_of_records);
            next_particles->flat_record_writer.reset(current.factor_variables->duplicate());
            next_particles->flat_record_writer->set_particle(NULL);
          }
          next_particles->flat_points.prepare(current.parameters, N*m);
        });
        size_t row = i*m + stored;
        if (number_of_records > 0)
        {
          records.resize(number_of_records);
          current.factor_variables->compact_values(records.data());
          next_particles->flat_points.set_point(row, current.parameters, records.data());
        }
        else
        {
          next_particles->flat_points.set_point(row, current.parameters);
        }
        ++stored;
      }

      if (tuned)
      {
        for (auto &t : this->tuners)
          t.second->score(i, before, current.parameters, this->particle_tunings[i]);
        current.tuning = NULL;
      }
      // as move(): the target over the SMC's index, which becomes the
      // particle's previous target (set_previous_target_evaluated_to_target_evaluated)
      double target = current.get_cached_likelihoods(this->index.get());
      const FactorVariables* factor_variables = current.factor_variables.get();
      int size = (factor_variables == NULL) ? -1 : factor_variables->adaptive_summary_size(this->index.get());
      std::call_once(final_setup, [&]()
      {
        next_particles->flat_final_points.prepare(current.parameters, N);
        if (size < 0)
          Rcpp::stop("SMCMCMCMove::move_flat: a factor has no adaptive summary.");
        summary_size = size;
        next_particles->adaptive_summaries.set_size(size_t(summary_size), N);
        next_particles->flat_summary_evaluator.reset(factor_variables->duplicate());
        next_particles->flat_summary_evaluator->set_particle(NULL);
        next_particles->adaptive_summary_evaluator = next_particles->flat_summary_evaluator.get();
        next_particles->adaptive_summary_index = this->index.get();
        next_particles->adaptive_summary_state = 1;
      });
      next_particles->flat_final_points.set_point(i, current.parameters);
      if (size != summary_size)
      {
        Rcpp::stop("SMCMCMCMove::move_flat: the particles' adaptive summaries differ in size.");
      }
      if (summary_size > 0)
        factor_variables->adaptive_summary(this->index.get(), next_particles->adaptive_summaries.colptr(i));
      next_particles->adaptive_previous_targets[i] = target;
    }
  };

  // the moves read the kernels and factors only (C98; checked by
  // supports_flat_moves, which admits exact and annealed-exact factors only)
  bool in_parallel = this->parallel && (this->factors != NULL) && this->factors->supports_parallel_chunks();
  if (this->parallel && !in_parallel)
    this->parallel_refused = true;
  run_chunk_groups(0, K, in_parallel, this->number_of_threads,
                   [&](size_t group, size_t first_chunk, size_t end_chunk)
                   {
                     Particle current;
                     Particle scratch;
                     std::vector<double> records;
                     for (size_t k = first_chunk; k < end_chunk; ++k)
                       run_chunk(k, current, scratch, records);
                   },
                   [&](size_t chunks_done)
                   {
                     if (this->progress_counts_iterations)
                       ProgressReporter::instance().update(double(std::min(N, chunks_done*B))*double(iterations));
                   });

  // as SMCWorker::move
  worker->set_seed(worker->get_seed() + N);
  this->rng->seed(worker->get_seed());
  next_particles->simulate_resampling_variables(*this->rng);
  worker->set_seed(worker->get_seed() + N);
}

void SMCMCMCMove::simulate_smc(SMCOutput* current_state)
{
  if (current_state->number_of_smc_iterations()==0)
  {
    // Simulate from the proposal: in chunks, kept flat with each point's
    // adaptive summary, when the moves that follow can start from a flat
    // generation (dev/plan.md 2.6a; the points are kept to be moved)
    if (this->can_simulate_flat())
    {
      this->simulate_generation_in_chunks(current_state, NULL, true, this->index.get(), true);
    }
    else
    {
      if (this->parallel)
        this->parallel_refused = true;
      this->simulate_proposal(current_state);
    }
  }
  else
  {
    // (dev/plan.md 3.21 (c2)) the factors' inner states that follow this
    // SMC's target take the random part of their step (rare-event
    // ABC-SMC^2), before resampling, as in arXiv:2211.02172's algorithm 4
    this->advance_inner_states(current_state);

    // update MCMC proposals
    this->mcmc->smc_adapt(current_state);
    
    current_state->normalise_and_resample_weights();
    // (c2') no move at a step without resampling, if asked
    if (this->move_only_after_resampling && !current_state->back().resampled_flag)
    {
      this->carry_generation(current_state);
      current_state->increment_smc_iteration();
      return;
    }
    // each particle's own tuning values (esjd, C211)
    this->prepare_tuning(current_state);
    
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Particles* current_particles = &current_state->back();
    Particles* next_particles = current_state->add_particles(current_particles);
    
    // (dev/plan.md 3.21 (c5)) the chains' iterations need not be stored when
    // the output is not written: only their final states are used (MCMC runs
    // and waste-free SMC use them all)
    this->mcmc->set_store_iterations(this->chains_only || (this->waste_free_chains > 0) || this->outputs_written);
    // (dev/plan.md 3.21 (c3)) adaptive_iterations: a pilot pass of one
    // iteration first
    this->start_adaptive_moves();
    if (this->waste_free_chains > 0)
    {
      if (this->adaptive_termination() != NULL)
        Rcpp::stop("SMC_with_MCMC: waste_free_chains needs a fixed number of MCMC iterations, not adaptive_iterations.");
      if (!this->can_move_flat())
        Rcpp::stop("SMC_with_MCMC: waste_free_chains needs exact or annealed factors, a move that does not adapt and a fixed number of MCMC iterations.");
      this->move_waste_free(next_particles,
                            current_particles);
    }
    else if (this->can_move_flat())
    {
      this->move_flat(next_particles,
                      current_particles);
    }
    else if (this->rich_chunk_mode() > 0)
    {
      this->move_rich_in_chunks(next_particles,
                                current_particles,
                                this->rich_chunk_mode() == 2);
    }
    else
    {
      if (this->parallel)
        this->parallel_refused = true;
      this->the_worker->move(next_particles,
                             current_particles);
    }
    this->finish_adaptive_moves(next_particles);
    
    current_state->increment_smc_iteration();
    
    // involves complete evaluation of weights using current adaptive param
  }
  
}

void SMCMCMCMove::evaluate_smc(SMCOutput* current_state)
{
  this->evaluate_smcfixed_part_smc(current_state);
  this->evaluate_smcadaptive_part_given_smcfixed_smc(current_state);
}

void SMCMCMCMove::evaluate_smcfixed_part_smc(SMCOutput* current_state)
{
  // a flat generation's smcfixed parts were evaluated when it was simulated
  if (current_state->back().is_flat())
    return;
  this->the_worker->smcfixed_weight(this->index.get(),
                                    current_state->back());
}

void SMCMCMCMove::evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* current_state)
{
  // set sequencer to have values from conditioned_on_parameters
  if (!this->sequencer_limit_is_fixed)
    Rcpp::stop("SMCMCMCMove::evaluate_smcadaptive_part_given_smcfixed_smc - need fixed sequencer limit if we are not conditioning on parameters.");
  
  // iterate until stop.
  bool terminate = FALSE;
  while (terminate==FALSE)
  {
    check_interrupt(); // between SMC steps (dev/plan.md 3.11)
    this->sequencer.find_desired_criterion(current_state);
    // (involves evaluating adaptive weights, using Sequencer)
    this->sequencer.find_next_target_bisection(current_state,this->index.get());
    this->finish_weighting_from_summaries(this->index.get(), current_state->back());
    if (this->show_progress && !this->progress_counts_iterations)
    {
      std::string status;
      double fraction = this->sequencer.progress_fraction(status);
      std::ostringstream text;
      text << "target " << current_state->number_of_smc_iterations() + 1 << ", " << status
           << ", ESS " << int(current_state->back().ess);
      ProgressReporter::instance().update((fraction < 0.0) ? 0.0 : 100.0*fraction, text.str());
    }
    
    current_state->log_likelihood = current_state->log_likelihood + current_state->calculate_latest_log_normalising_constant_ratio();
    // One log-likelihood and one time per stored particle set: recorded here,
    // before the termination check, so that the particle set written out by
    // the final MCMC move below has both.
    current_state->set_llhd(current_state->log_likelihood);
    current_state->set_time_and_reset_start();
    
    current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();

    // check termination, using sequencer
    if (this->sequencer.check_termination())
    {
      
      if (this->mcmc_at_last_step)
      {
        this->simulate_smc(current_state);
        current_state->set_llhd(current_state->log_likelihood);
        current_state->set_time_and_reset_start();
      }
      current_state->terminate();
      break;
    }
    
    this->simulate_smc(current_state);
    current_state->back().set_previous_target_evaluated_to_target_evaluated();
  }
}

MoveOutput* SMCMCMCMove::move(RandomNumberGenerator &rng,
                              const Particle &particle) const
{
  // The MCMC moves read their acceptance-ratio denominators from the
  // particle's per-factor cache, so every factor in the target must have been
  // evaluated at the current parameters. The SMC has evaluated the factors in
  // its own index; evaluate the ones it skipped (cancelled_index) here. Only
  // those are evaluated, so random (pseudo-marginal) factors the SMC has
  // already estimated are not re-estimated.
  // A chain of unadjusted moves only never reads cached values, and the user
  // may not have supplied the functions to evaluate the factors, so nothing is
  // evaluated here (the behaviour before the cache was used).
  if (!this->mcmc->reads_cached_factors())
    return this->mcmc->run(rng, particle);
  
  MoveOutput* output;
  if (this->cancelled_index!=NULL)
  {
    Particle start = particle;
    start.evaluate_likelihoods(this->cancelled_index.get());
    output = this->mcmc->run(rng, start);
  }
  else
  {
    output = this->mcmc->run(rng, particle);
  }
  
  // The SMC uses target_evaluated (over its own index) in its next weight
  // update, so set it from the cache rather than leaving whatever the last
  // move happened to store there.
  output->back().target_evaluated = output->back().get_cached_likelihoods(this->index.get());
  return output;
}

SMCOutput* SMCMCMCMove::specific_run(const Parameters &conditioned_on_parameters)
{
  SMCOutput* simulation = this->initialise_smc(conditioned_on_parameters);
  this->simulate_smc(simulation, conditioned_on_parameters);
  this->evaluate_smc(simulation, conditioned_on_parameters);
  simulation->normalise_and_resample_weights();
  return simulation;
}

SMCOutput* SMCMCMCMove::specific_initialise_smc(const Parameters &conditioned_on_parameters)
{
  SMCOutput* output = new SMCOutput(this, this->lag, this->lag_proposed, this->results_name);
  return output;
}

void SMCMCMCMove::simulate_smc(SMCOutput* current_state,
                               const Parameters &conditioned_on_parameters)
{
  if (current_state->number_of_smc_iterations()==0)
  {
    // Simulate from the proposal.
    this->simulate_proposal(current_state, conditioned_on_parameters);
  }
  else
  {
    this->advance_inner_states(current_state); // (dev/plan.md 3.21 (c2))

    // update MCMC proposals
    this->mcmc->smc_adapt(current_state);
    
    current_state->normalise_and_resample_weights();
    // (c2') no move at a step without resampling, if asked
    if (this->move_only_after_resampling && !current_state->back().resampled_flag)
    {
      this->carry_generation(current_state);
      current_state->increment_smc_iteration();
      return;
    }
    
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Particles* current_particles = &current_state->back();
    Particles* next_particles = current_state->add_particles(current_particles);
    this->mcmc->set_store_iterations(this->chains_only || this->outputs_written); // (c5: stored only when written)
    this->start_adaptive_moves(); // (dev/plan.md 3.21 (c3))
    this->the_worker->move(next_particles,
                           current_particles);
    this->finish_adaptive_moves(next_particles);
    
    current_state->increment_smc_iteration();
    // involves complete evaluation of weights using current adaptive param
  }
  
}

void SMCMCMCMove::evaluate_smc(SMCOutput* current_state,
                               const Parameters &conditioned_on_parameters)
{
  this->evaluate_smcfixed_part_smc(current_state,
                                   conditioned_on_parameters);
  this->evaluate_smcadaptive_part_given_smcfixed_smc(current_state,
                                                     conditioned_on_parameters);
}

void SMCMCMCMove::evaluate_smcfixed_part_smc(SMCOutput* current_state,
                                             const Parameters &conditioned_on_parameters)
{
  
  this->the_worker->smcfixed_weight(this->index.get(),
                                    current_state->back());
}

void SMCMCMCMove::evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* current_state,
                                                               const Parameters &conditioned_on_parameters)
{
  if (this->follows_caller_target)
  {
    // (dev/plan.md 3.21 (c2)) a new output is run through the caller's
    // targets so far; afterwards an evaluation is a trial at the caller's
    // current (candidate) target: the log likelihood at the last target it
    // was advanced to, times the estimate of the ratio from its weighted
    // points (no random numbers, so the caller's bisection may try any target)
    if (!current_state->caught_up)
    {
      this->follow_from_start(current_state, conditioned_on_parameters);
      return;
    }
    // (an estimate of zero stays zero: its points cannot be reweighted)
    double ratio = (current_state->committed_log_likelihood == -arma::datum::inf) ? 0.0 :
      this->step_to_target(current_state, this->followed_target(conditioned_on_parameters));
    current_state->log_likelihood = current_state->committed_log_likelihood + ratio;
    current_state->set_latest_llhd(current_state->log_likelihood); // (the latest generation's, as weighted now)
    return;
  }
  // set sequencer to have values from conditioned_on_parameters
  if (!this->sequencer_limit_is_fixed)
    this->sequencer.set_next_with_parameter(conditioned_on_parameters);
  else
    this->sequencer.reset();
  
  // iterate until stop.
  bool terminate = FALSE;
  while (terminate==FALSE)
  {
    check_interrupt(); // between SMC steps (dev/plan.md 3.11)
    this->sequencer.find_desired_criterion(current_state);
    
    // (involves evaluating adaptive weights, using Sequencer)
    this->sequencer.find_next_target_bisection(current_state,
                                               this->index.get());
    this->finish_weighting_from_summaries(this->index.get(), current_state->back());
    current_state->log_likelihood = current_state->log_likelihood + current_state->calculate_latest_log_normalising_constant_ratio();
    current_state->set_llhd(current_state->log_likelihood);
    
    current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();
    
    current_state->set_time_and_reset_start();

    // (dev/plan.md 3.21 (c4)) early stopping in an MCMC move's proposal
    if (!std::isnan(this->run_threshold) &&
        (current_state->log_likelihood + this->remaining_maximum(conditioned_on_parameters) < this->run_threshold))
    {
      this->run_stopped = true;
      current_state->terminate();
      break;
    }

    // check termination, using sequencer
    if (this->sequencer.check_termination())
    {
      if (this->mcmc_at_last_step)
      {
        this->simulate_smc(current_state,conditioned_on_parameters);
        current_state->set_llhd(current_state->log_likelihood);
        current_state->set_time_and_reset_start();
      }
      current_state->terminate();
      break;
    }
    
    this->simulate_smc(current_state, conditioned_on_parameters);
    
    current_state->back().set_previous_target_evaluated_to_target_evaluated();
  }
}

void SMCMCMCMove::subsample_simulate_smc(SMCOutput* current_state)
{
  if (current_state->number_of_smc_iterations()==0)
  {
    // Simulate from the proposal.
    this->simulate_proposal(current_state);
  }
  else
  {
    // update MCMC proposals
    this->mcmc->smc_adapt(current_state);
    
    current_state->normalise_and_resample_weights();
    
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Particles* current_particles = &current_state->back();
    Particles* next_particles = current_state->add_particles(current_particles);
    this->the_worker->subsample_move(next_particles,
                                     current_particles);
    
    current_state->increment_smc_iteration();
    // involves complete evaluation of weights using current adaptive param
  }
  
}

void SMCMCMCMove::subsample_simulate_smc(SMCOutput* current_state,
                                         const Parameters &conditioned_on_parameters)
{
  if (current_state->number_of_smc_iterations()==0)
  {
    // Simulate from the proposal.
    this->simulate_proposal(current_state, conditioned_on_parameters);
  }
  else
  {
    // update MCMC proposals
    this->mcmc->smc_adapt(current_state);
    
    current_state->normalise_and_resample_weights();
    
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Particles* current_particles = &current_state->back();
    Particles* next_particles = current_state->add_particles(current_particles);
    this->the_worker->subsample_move(next_particles,
                                     current_particles);
    
    current_state->increment_smc_iteration();
    // involves complete evaluation of weights using current adaptive param
  }
  
}

void SMCMCMCMove::subsample_evaluate_smc(SMCOutput* current_state)
{
  this->subsample_evaluate_smcfixed_part_smc(current_state);
  this->subsample_evaluate_smcadaptive_part_given_smcfixed_smc(current_state);
}

void SMCMCMCMove::subsample_evaluate_smc(SMCOutput* current_state,
                                         const Parameters &conditioned_on_parameters)
{
  this->subsample_evaluate_smcfixed_part_smc(current_state,
                                             conditioned_on_parameters);
  this->subsample_evaluate_smcadaptive_part_given_smcfixed_smc(current_state,
                                                               conditioned_on_parameters);
}

void SMCMCMCMove::subsample_evaluate_smcfixed_part_smc(SMCOutput* current_state)
{
  this->the_worker->subsample_smcfixed_weight(this->index.get(),
                                              current_state->back());
}

void SMCMCMCMove::subsample_evaluate_smcfixed_part_smc(SMCOutput* current_state,
                                                       const Parameters &conditioned_on_parameters)
{
  this->the_worker->subsample_smcfixed_weight(this->index.get(),
                                              current_state->back());
}

void SMCMCMCMove::subsample_evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* current_state)
{
  // set sequencer to have values from conditioned_on_parameters
  
  // need to do something with sequencer for each particle
  
  // either set sequencer with paraneters - for now, just setting an end point is fine, but could make it more sophisticated
  
  // or reset sequencer
  
  // iterate until stop.
  bool terminate = FALSE;
  while (terminate==FALSE)
  {
    check_interrupt(); // between SMC steps (dev/plan.md 3.11)
    this->sequencer.find_desired_criterion(current_state);
    
    // (involves evaluating adaptive weights, using Sequencer)
    this->sequencer.subsample_find_next_target_bisection(current_state,
                                                         this->index.get());
    
    current_state->subsample_log_likelihood = current_state->subsample_log_likelihood + current_state->calculate_latest_log_normalising_constant_ratio();
    
    current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();
    
    current_state->set_llhd(current_state->subsample_log_likelihood);
    current_state->set_time_and_reset_start();

    // check termination, using sequencer
    if (this->sequencer.check_termination())
    {
      if (this->mcmc_at_last_step)
      {
        // The final move targets the same (subsampled) distribution as the
        // other moves of this run.
        this->subsample_simulate_smc(current_state);
        current_state->set_llhd(current_state->subsample_log_likelihood);
        current_state->set_time_and_reset_start();
      }
      current_state->terminate();
      break;
    }
    
    this->subsample_simulate_smc(current_state);
    
    current_state->back().subsample_set_previous_target_evaluated_to_target_evaluated();
  }
}

void SMCMCMCMove::subsample_evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* current_state,
                                                                         const Parameters &conditioned_on_parameters)
{
  // set sequencer to have values from conditioned_on_parameters
  if (!this->sequencer_limit_is_fixed)
    this->sequencer.set_next_with_parameter(conditioned_on_parameters);
  else
    this->sequencer.reset();
  
  // need to do something with sequencer for each particle
  
  // either set sequencer with paraneters - for now, just setting an end point is fine, but could make it more sophisticated
  
  // or reset sequencer
  
  // iterate until stop.
  bool terminate = FALSE;
  while (terminate==FALSE)
  {
    check_interrupt(); // between SMC steps (dev/plan.md 3.11)
    this->sequencer.find_desired_criterion(current_state);
    
    // (involves evaluating adaptive weights, using Sequencer)
    this->sequencer.subsample_find_next_target_bisection(current_state,
                                                         this->index.get());
    
    current_state->subsample_log_likelihood = current_state->subsample_log_likelihood + current_state->calculate_latest_log_normalising_constant_ratio();
    
    current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();
    
    current_state->set_llhd(current_state->subsample_log_likelihood);
    current_state->set_time_and_reset_start();

    // check termination, using sequencer
    if (this->sequencer.check_termination())
    {
      if (this->mcmc_at_last_step)
      {
        // The final move targets the same (subsampled) distribution as the
        // other moves of this run.
        this->subsample_simulate_smc(current_state, conditioned_on_parameters);
        current_state->set_llhd(current_state->subsample_log_likelihood);
        current_state->set_time_and_reset_start();
      }
      current_state->terminate();
      break;
    }
    
    this->subsample_simulate_smc(current_state,
                                 conditioned_on_parameters);
    
    current_state->back().subsample_set_previous_target_evaluated_to_target_evaluated();
  }
}

void SMCMCMCMove::weight_for_adapting_sequence(const Index* index,
                                               Particles &current_particles)
{
  // from the generation's adaptive summaries when available (dev/plan.md
  // 2.6a, G2; the bisection calls this once per trial target)
  this->weight_from_summaries_if_possible(this->index.get(),
                                          current_particles);
}

MoveOutput* SMCMCMCMove::subsample_move(RandomNumberGenerator &rng,
                                        const Particle &particle) const
{
  return this->mcmc->subsample_run(rng,
                                   particle);
}

void SMCMCMCMove::subsample_weight_for_adapting_sequence(const Index* index,
                                                         Particles &current_particles)
{
  this->the_worker->subsample_smcadaptive_given_smcfixed_weight(this->index.get(),
                                                                current_particles);
}


// --- rare-event ABC-SMC^2: the inner SMC following the caller's target
// (dev/plan.md 3.21 (c2); arXiv:2211.02172, algorithms 2 and 4) -------------

void SMCMCMCMove::set_follows_caller_target(const std::string &variable)
{
  this->follows_caller_target = true;
  this->followed_variable = variable;
  // the caller re-evaluates this factor at each of its targets
  this->smcfixed_flag = false;
  // each point is reweighted and moved through its own particle
  this->force_rich_generations = true;
}

double SMCMCMCMove::followed_target(const Parameters &conditioned_on_parameters) const
{
  return conditioned_on_parameters[this->followed_variable][0];
}

// reweight the latest generation from its previous target to `target` (the
// points read it through the sequencer's shared schedule parameters, merged
// after the caller's parameters): the weights start again from the step's
// starting weights on each call, so trials leave nothing behind; returns the
// log of the estimated ratio of likelihoods
double SMCMCMCMove::step_to_target(SMCOutput *current_state, double target)
{
  this->sequencer.schedule_parameters[this->followed_variable] = target;
  this->weight_for_adapting_sequence(this->index.get(), current_state->back());
  current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
  return current_state->calculate_latest_log_normalising_constant_ratio();
}

// the entries of a generation just made (one log likelihood and one time per
// generation, as SMCOutput's writing expects); a later reweighting of the same
// generation replaces its log likelihood
static void record_generation(SMCOutput *current_state)
{
  current_state->set_llhd(current_state->committed_log_likelihood);
  current_state->set_time_and_reset_start();
}

// a step at the current target: reweight to it (the trial the caller has just
// made, now kept), then resample and move u (the MCMC kernel at this target)
void SMCMCMCMove::advance_at_current_target(SMCOutput *current_state, const Parameters &conditioned_on_parameters)
{
  double target = this->followed_target(conditioned_on_parameters);
  // an estimate of zero (every point outside the kernel) stays zero, and its
  // weights cannot be normalised: the state stops there (the caller's
  // particle has weight zero)
  if (current_state->committed_log_likelihood != -arma::datum::inf)
  {
    current_state->committed_log_likelihood += this->step_to_target(current_state, target);
    current_state->back().set_previous_target_evaluated_to_target_evaluated();
    current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();
    current_state->set_latest_llhd(current_state->committed_log_likelihood);
    if (current_state->committed_log_likelihood != -arma::datum::inf)
    {
      this->simulate_smc(current_state, conditioned_on_parameters);
      current_state->back().set_previous_target_evaluated_to_target_evaluated();
      record_generation(current_state);
    }
  }
  current_state->log_likelihood = current_state->committed_log_likelihood;
}

// the caller's move stage: the caller's chosen target joins the targets so
// far (every caller particle is at the same one), and this output steps to it
void SMCMCMCMove::advance_at_caller_target_smc(SMCOutput *current_state, const Parameters &conditioned_on_parameters)
{
  if (!this->follows_caller_target || !current_state->caught_up)
    return;
  double target = this->followed_target(conditioned_on_parameters);
  if (this->followed_targets.empty() || (this->followed_targets.back() != target))
    this->followed_targets.push_back(target);
  this->advance_at_current_target(current_state, conditioned_on_parameters);
}

// a new output (its points just drawn from u's distribution given theta): the
// same steps as an output that has followed the caller from its start, so
// that both have the same law (the pseudo-marginal ratio of the caller's
// moves needs this): reweighted, resampled and moved at each of the targets
// the caller has advanced to so far (none, for an output made at the
// caller's start: it is then the draws themselves, unweighted). Its value is
// then a trial at the caller's current target if it is not among those —
// as for an output made at the caller's start, which the caller evaluates at
// a first target its sequencer may still change (an adaptive schedule
// searches from there in either direction), so nothing is kept at it.
void SMCMCMCMove::follow_from_start(SMCOutput *current_state, const Parameters &conditioned_on_parameters)
{
  double current = this->followed_target(conditioned_on_parameters);
  current_state->committed_log_likelihood = 0.0;
  record_generation(current_state); // (the first generation's entries)
  for (size_t k = 0; k < this->followed_targets.size(); ++k)
  {
    current_state->committed_log_likelihood += this->step_to_target(current_state, this->followed_targets[k]);
    current_state->back().set_previous_target_evaluated_to_target_evaluated();
    current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();
    current_state->set_latest_llhd(current_state->committed_log_likelihood);
    if (current_state->committed_log_likelihood == -arma::datum::inf)
      break; // (an estimate of zero stays zero)
    this->simulate_smc(current_state, conditioned_on_parameters);
    current_state->back().set_previous_target_evaluated_to_target_evaluated();
    record_generation(current_state);
  }
  current_state->caught_up = true;
  current_state->log_likelihood = current_state->committed_log_likelihood;
  if ((current_state->committed_log_likelihood != -arma::datum::inf) &&
      (this->followed_targets.empty() || (this->followed_targets.back() != current)))
  {
    current_state->log_likelihood += this->step_to_target(current_state, current);
    current_state->set_latest_llhd(current_state->log_likelihood);
  }
}

// at this SMC's move stage: the inner states of the particles' factors
void SMCMCMCMove::advance_inner_states(SMCOutput *current_state)
{
  if ((this->factors == NULL) || !this->factors->advances_at_caller_target())
    return;
  Particles &generation = current_state->back();
  for (size_t i = 0; i < generation.size(); ++i)
  {
    Particle &particle = generation[i]->back();
    if (particle.factor_variables != NULL)
      particle.factor_variables->advance_at_caller_target();
  }
}


// --- adaptive number of MCMC iterations (dev/plan.md 3.21 (c3); South,
// Pettitt & Drovandi 2019, as in arXiv:2211.02172 §3.3.3) --------------------

AdaptiveIterationsMCMCTermination* SMCMCMCMove::adaptive_termination() const
{
  if (this->mcmc == NULL)
    return NULL;
  return dynamic_cast<AdaptiveIterationsMCMCTermination*>(this->mcmc->get_termination());
}

void SMCMCMCMove::start_adaptive_moves()
{
  AdaptiveIterationsMCMCTermination* termination = this->adaptive_termination();
  if (termination == NULL)
    return;
  if (this->chains_only)
    Rcpp::stop("adaptive_iterations chooses the number of iterations of each SMC step from the particles: "
               "it is for SMC moves, not MCMC (use iterations(n)).");
  termination->set_iterations(1);
}

// TRUE if a particle's own (non-fixed) variables differ from its ancestor's
static bool particle_moved(const Particle &now, const Particle &before)
{
  for (const std::string &name : now.parameters.get_nonfixed_vector_variables())
  {
    const arma::mat &x = now.parameters[name];
    const arma::mat &y = before.parameters[name];
    if ((x.n_rows != y.n_rows) || (x.n_cols != y.n_cols) || arma::any(arma::vectorise(x != y)))
      return true;
  }
  return false;
}

// after the pilot pass: the acceptance rate p (the fraction of particles that
// moved), then R - 1 more iterations from each particle, R = ceil(log c /
// log(1 - p)) (sequentially, with their own seeds; each particle's output is
// replaced by the longer chain's, whose start was the pilot's state)
void SMCMCMCMove::finish_adaptive_moves(Particles* next_particles)
{
  AdaptiveIterationsMCMCTermination* termination = this->adaptive_termination();
  if (termination == NULL)
    return;
  size_t N = next_particles->size();
  size_t moved = 0;
  for (size_t i = 0; i < N; ++i)
  {
    const Particle &now = (*next_particles)[i]->back();
    if ((now.previous_self != NULL) && particle_moved(now, *now.previous_self))
      ++moved;
  }
  double acceptance_rate = (N == 0) ? 0.0 : double(moved)/double(N);
  size_t R = termination->iterations_for(acceptance_rate);
  if (R > 1)
  {
    termination->set_iterations(R - 1);
    RandomNumberGenerator local_rng(*this->rng);
    size_t seed = this->the_worker->get_seed();
    local_rng.seed(seed);
    for (size_t i = 0; i < N; ++i)
    {
      Particle &start = (*next_particles)[i]->back();
      const Particle* ancestor = start.previous_self;
      MoveOutput* longer = this->move(local_rng, start);
      longer->back().previous_self = ancestor;
      next_particles->particles[i].reset(longer);
    }
    this->the_worker->set_seed(seed + N);
  }
  termination->set_iterations(1);
}


// --- early stopping in rare-event ABC-MCMC (dev/plan.md 3.21 (c4)) ----------

// only for ABC kernels (the lowering's flag) and tolerances that never increase
void SMCMCMCMove::set_stops_early(bool kernels_only)
{
  this->stops_early = false;
  if (!kernels_only || this->sequencer.schedules.empty())
    return;
  const std::vector<double> &schedule = this->sequencer.schedules.back();
  for (size_t k = 1; k < schedule.size(); ++k)
    if (schedule[k] > schedule[k - 1])
      return;
  this->stops_early = true;
}

// log max K_T - log max K_s: the factors' maxima at the final tolerance and at
// the current one (+Inf if one is not known: no stop)
double SMCMCMCMove::remaining_maximum(const Parameters &conditioned_on_parameters) const
{
  if ((this->factors == NULL) || this->sequencer.variable_names.empty() || this->sequencer.schedules.empty())
    return arma::datum::inf;
  const std::string &variable = this->sequencer.variable_names.back();
  // (the schedule's value first: merge keeps the first set's value of a name)
  Parameters at_current = this->sequencer.schedule_parameters.deep_copy().merge(conditioned_on_parameters);
  double current = this->factors->log_maximum(this->index.get(), at_current);
  Parameters at_final = at_current.deep_copy();
  at_final[variable] = this->sequencer.schedules.back().back();
  double final_maximum = this->factors->log_maximum(this->index.get(), at_final);
  if (!std::isfinite(current) || !std::isfinite(final_maximum))
    return arma::datum::inf;
  return final_maximum - current;
}

double SMCMCMCMove::run_with_threshold(SMCOutput* simulation,
                                       const Parameters &conditioned_on_parameters,
                                       double threshold,
                                       bool &below)
{
  this->run_threshold = threshold;
  this->run_stopped = false;
  this->simulate_smc(simulation, conditioned_on_parameters);
  double value = simulation->evaluate(conditioned_on_parameters);
  this->run_threshold = arma::datum::nan;
  below = this->run_stopped;
  this->run_stopped = false;
  return below ? -arma::datum::inf : value;
}


// (dev/plan.md 3.21 (c2')) the next generation without moves: each particle
// copied (its weight carried on; not resampled, so its own ancestor), and
// the next resampling's random numbers drawn as a move would (the same seeds)
void SMCMCMCMove::carry_generation(SMCOutput* current_state)
{
  Particles* current_particles = &current_state->back();
  Particles* next_particles = current_state->add_particles(current_particles);
  size_t N = current_particles->size();
  for (size_t i = 0; i < N; ++i)
  {
    const Particle &ancestor = (*current_particles)[current_particles->ancestor_variables[i]]->back();
    next_particles->push_back(new SinglePointMoveOutput(ancestor));
    next_particles->back()->back().previous_self = &ancestor;
  }
  this->the_worker->set_seed(this->the_worker->get_seed() + this->the_worker->get_number_of_particles());
  this->rng->seed(this->the_worker->get_seed());
  next_particles->simulate_resampling_variables(*this->rng);
  this->the_worker->set_seed(this->the_worker->get_seed() + this->the_worker->get_number_of_particles());
}

}
