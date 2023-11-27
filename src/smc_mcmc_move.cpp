#include <numeric>
#include "smc_mcmc_move.h"
#include "smc_output.h"
#include "cess_smc_criterion.h"
#include "exact_likelihood_estimator.h"
#include "annealed_likelihood_estimator.h"
#include "utils.h"
#include "parameter_particle_simulator.h"
//#include "rcppparallel_smc_worker.h"
#include "sequential_smc_worker.h"
#include "mcmc.h"
#include "move_output.h"
#include "standard_mcmc_output.h"
#include "custom_distribution_proposal_kernel.h"
#include "vector_factors.h"
#include "vector_single_index.h"

SMCMCMCMove::SMCMCMCMove()
   :SMC()
{
  this->mcmc = NULL;
  this->index = NULL;
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
     std::max<size_t>(2,lag_in),
     lag_proposed_in,
     mcmc_in->get_proposals(),
     -1.0,
     false,
     true,
     true,
     transform_proposed_particles,
     results_name_in)
{
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
  this->index = new VectorSingleIndex(indices);

  this->factors = new VectorFactors(likelihood_estimators);
  
  this->particle_simulator = NULL;
  
  this->initial_particles = initial_points_in;
  this->log_probabilities_of_initial_values = log_probabilities_of_initial_values_in;
  this->proposed_particles_inputted = true;
  
  if (initial_points_in.size()==0)
    Rcpp::stop("SMCMCMCMove::SMCMCMCMove - need more than one initialising point.");
  
  //this->model_and_algorithm.particle_simulator = new ParameterParticleSimulator(simulate_proposal_in,
  //                                                                              this->model_and_algorithm.likelihood_estimators);
  
  if (parallel_in==TRUE)
  {
  }
  else
  {
    this->the_worker = new SequentialSMCWorker(this);
  }
  
  std::vector<double> schedule_in;
  schedule_in.push_back(0.0);
  schedule_in.push_back(1.0);
  std::string variable_in = "power";
  SMCCriterion* smc_criterion = new CESSSMCCriterion(-1.0);
  this->sequencer = Sequencer(this->the_worker,
                              schedule_in,
                              variable_in,
                              25,
                              smc_criterion);
  
  this->mcmc = mcmc_in;
  this->mcmc->set_index(new VectorSingleIndex(indices));
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
     std::max<size_t>(2,lag_in),
     lag_proposed_in,
     mcmc_in->get_proposals(),
     -1.0,
     true,
     true,
     true,
     transform_proposed_particles,
     results_name_in)
{
  mcmc_in->set_proposal_parameters(&this->algorithm_parameters);
  proposal_in->set_proposal_parameters(&this->algorithm_parameters);
  
  //std::vector<size_t> indices;
  //indices.reserve(likelihood_estimators_in.size());
  //for (size_t i=0; i<likelihood_estimators_in.size(); ++i)
  //  indices.push_back(i);
  //this->index = new VectorSingleIndex(indices);
  this->index = without_cancelled_index;
  
  this->factors = new VectorFactors(likelihood_estimators_in);
  
  this->particle_simulator = new ParameterParticleSimulator(proposal_in,
                                                            likelihood_estimators_in);
  
  this->proposed_particles_inputted = false;

  //this->model_and_algorithm.particle_simulator = new ParameterParticleSimulator(simulate_proposal_in,
  //                                                                              this->model_and_algorithm.likelihood_estimators);
  
  if (parallel_in==TRUE)
  {
    //this->the_worker = new RcppParallelSMCWorker(this,grain_size_in);
    Rcpp::stop("Parallel worker not set up.");
  }
  else
  {
    this->the_worker = new SequentialSMCWorker(this);
  }
  
  std::vector<double> schedule_in;
  schedule_in.push_back(0.0);
  schedule_in.push_back(1.0);
  std::string variable_in = "power";
  SMCCriterion* smc_criterion = new CESSSMCCriterion(-1.0);
  this->sequencer = Sequencer(this->the_worker,
                              schedule_in,
                              variable_in,
                              25,
                              smc_criterion);
  
  this->mcmc = mcmc_in;
  if (this->mcmc!=NULL)
    this->mcmc->set_index_if_null(full_index);
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
     std::max<size_t>(2,lag_in),
     lag_proposed_in,
     mcmc_in->get_proposals(),
     -1.0,
     false,
     true,
     true,
     transform_proposed_particles,
     results_name_in)
{
  mcmc_in->set_proposal_parameters(&this->algorithm_parameters);
  
  //std::vector<size_t> indices;
  //indices.reserve(likelihood_estimators_in.size());
  //for (size_t i=0; i<likelihood_estimators_in.size(); ++i)
  //  indices.push_back(i);
  //this->index = new VectorSingleIndex(indices);
  this->index = without_cancelled_index;
  
  this->factors = new VectorFactors(likelihood_estimators_in);
  
  this->particle_simulator = NULL;
  
  this->initial_particles = initial_points_in;
  
  //arma::mat thing = this->initial_particles[0]["tau"];
  this->log_probabilities_of_initial_values = log_probabilities_of_initial_values_in;
  this->proposed_particles_inputted = true;
  
  if (initial_points_in.size()==0)
    Rcpp::stop("SMCMCMCMove::SMCMCMCMove - need more than one initialising point.");
  
  //this->model_and_algorithm.particle_simulator = new ParameterParticleSimulator(simulate_proposal_in,
  //                                                                              this->model_and_algorithm.likelihood_estimators);
  
  if (parallel_in==TRUE)
  {
    //this->the_worker = new RcppParallelSMCWorker(this,grain_size_in);
    Rcpp::stop("Parallel worker not set up.");
  }
  else
  {
    this->the_worker = new SequentialSMCWorker(this);
  }
  
  std::vector<double> schedule_in;
  schedule_in.push_back(0.0);
  schedule_in.push_back(1.0);
  std::string variable_in = "power";
  SMCCriterion* smc_criterion = new CESSSMCCriterion(-1.0);
  this->sequencer = Sequencer(this->the_worker,
                              schedule_in,
                              variable_in,
                              100,
                              smc_criterion);
  
  this->mcmc = mcmc_in;
  if (this->mcmc!=NULL)
    this->mcmc->set_index_if_null(full_index);
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
:SMC(rng_in, seed_in, data_in, algorithm_parameters, number_of_particles_in, std::max<size_t>(2,lag_in), lag_proposed_in, mcmc_in->get_proposals(), -1.0, true, true, true, transform_proposed_particles, results_name_in)
{
  mcmc_in->set_proposal_parameters(&this->algorithm_parameters);
  
  IndependentProposalKernel* proposal = new CustomDistributionProposalKernel(simulate_proposal_in,
                                                                             evaluate_log_proposal_in);
  proposal->set_proposal_parameters(&this->algorithm_parameters);
  //Parameters candidate_parameters = proposal->independent_simulate(*this->rng);
  
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
  this->index = new VectorSingleIndex(indices);
  
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
  
  if (parallel_in==TRUE)
  {
  }
  else
  {
    this->the_worker = new SequentialSMCWorker(this);
  }
  
  std::vector<double> schedule_in;
  schedule_in.push_back(0.0);
  schedule_in.push_back(1.0);
  std::string variable_in = "power";
  SMCCriterion* smc_criterion = new CESSSMCCriterion(-1.0);
  this->sequencer = Sequencer(this->the_worker,
                              schedule_in,
                              variable_in,
                              25,
                              smc_criterion);
  
  this->mcmc = mcmc_in;
  this->mcmc->set_index(new VectorSingleIndex(indices));
  this->mcmc_at_last_step = mcmc_at_last_step_in;
}

/*
void SMCMCMCMove::set_multiple_mcmc(RandomNumberGenerator* rng_in,
                                    size_t* seed_in,
                                    const Data* data_in,
                                    MCMC* mcmc_in,
                                    SimulateDistributionPtr simulate_proposal_in,
                                    EvaluateLogDistributionPtr evaluate_log_prior_in,
                                    EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
                                    bool parallel_in,
                                    size_t grain_size_in)
{
  this->model_and_algorithm.likelihood_estimators.resize(0);
  this->model_and_algorithm.likelihood_estimators.reserve(1);
  this->model_and_algorithm.likelihood_estimators.push_back(new ExactLikelihoodEstimator(rng_in,
                                                                                         seed_in,
                                                                                         data_in,
                                                                                         evaluate_log_prior_in,
                                                                                         evaluate_log_likelihood_in,
                                                                                         true));
  
  this->model_and_algorithm.particle_simulator = new ParameterParticleSimulator(simulate_proposal_in,
                                                                                this->model_and_algorithm.likelihood_estimators);
  
  if (parallel_in==TRUE)
  {
  }
  else
  {
    this->the_worker = new SequentialSMCWorker(this,
                                               this->model_and_algorithm.particle_simulator,
                                               this->model_and_algorithm.likelihood_estimators);
  }
  
  std::vector<double> schedule_in;
  schedule_in.push_back(0.0);
  schedule_in.push_back(1.0);
  std::string variable_in = "power";
  this->model_and_algorithm.smc_criterion = new CESSSMCCriterion(-1.0);
  this->model_and_algorithm.smc_termination = NULL;
  this->sequencer = Sequencer(this->the_worker,
                              schedule_in,
                              variable_in,
                              this->model_and_algorithm.smc_criterion,
                              this->model_and_algorithm.smc_termination);
  
  this->mcmc = mcmc_in;
}
*/

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
  :SMC(rng_in, seed_in, data_in, algorithm_parameters, number_of_particles_in, std::max<size_t>(2,lag_in), lag_proposed_in, mcmc_in->get_proposals(), resampling_desired_ess_in, true, true, true, transform_proposed_particles, results_name_in)
{
  mcmc_in->set_proposal_parameters(&this->algorithm_parameters);
  
  std::string variable_in = "power";
  
  // Need to construct LikelihoodEstimator to read in to this constructor.
  IndependentProposalKernel* proposal = new CustomDistributionProposalKernel(simulate_proposal_in,
                                                                             evaluate_log_proposal_in);
  proposal->set_proposal_parameters(&this->algorithm_parameters);
  //Parameters candidate_parameters = proposal->independent_simulate(*this->rng);
  
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
  this->index = new VectorSingleIndex(indices);
  
  /*
  for (auto i=likelihood_estimators.begin();
       i!=likelihood_estimators.end();
       ++i)
  {
    (*i)->setup(candidate_parameters);
  }
  */
  
  this->factors = new VectorFactors(likelihood_estimators);
  
  // Need to construct LikelihoodEstimator to read in to this constructor.
  this->particle_simulator = new ParameterParticleSimulator(proposal,
                                                            likelihood_estimators);

  if (parallel_in==true)
  {
  }
  else
  {
    this->the_worker = new SequentialSMCWorker(this);
  }
  
  std::vector<double> schedule_in;
  schedule_in.push_back(0.0);
  schedule_in.push_back(1.0);
  SMCCriterion* smc_criterion = new CESSSMCCriterion(annealing_desired_cess_in);
  //this->model_and_algorithm.smc_termination = NULL;
  this->sequencer = Sequencer(this->the_worker,
                              schedule_in,
                              variable_in,
                              number_of_bisections_in,
                              smc_criterion);
  //this->sequencer_parameters = &this->sequencer.schedule_parameters;
  
  this->mcmc = mcmc_in;
  this->mcmc->set_index(new VectorSingleIndex(indices));
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
:SMC(rng_in, seed_in, data_in, algorithm_parameters, number_of_particles_in, std::max<size_t>(2,lag_in), lag_proposed_in, mcmc_in->get_proposals(), resampling_desired_ess_in, true, true, true, transform_proposed_particles, results_name_in)
{
  mcmc_in->set_proposal_parameters(&this->algorithm_parameters);
  std::string variable_in = "power";
  // Need to construct LikelihoodEstimator to read in to this constructor.
  IndependentProposalKernel* proposal = new CustomDistributionProposalKernel(simulate_proposal_in,
                                                                             evaluate_log_proposal_in);
  proposal->set_proposal_parameters(&this->algorithm_parameters);
  //Parameters candidate_parameters = proposal->independent_simulate(*this->rng);
  
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
  this->index = new VectorSingleIndex(indices);
  
  this->factors = new VectorFactors(likelihood_estimators);
  
  /*
  for (auto i=likelihood_estimators.begin();
       i!=likelihood_estimators.end();
       ++i)
  {
    (*i)->setup(candidate_parameters);
  }
  */
  
  // Need to construct LikelihoodEstimator to read in to this constructor.
  this->particle_simulator = new ParameterParticleSimulator(proposal,
                                                            likelihood_estimators);
  
  if (parallel_in==true)
  {
  }
  else
  {
    this->the_worker = new SequentialSMCWorker(this);
  }
  
  std::vector<double> schedule_in;
  schedule_in.reserve(temperatures_in.size()+1);
  schedule_in.push_back(0.0);
  for (size_t i=0; i<temperatures_in.size(); ++i)
  {
    schedule_in.push_back(temperatures_in[i]);
  }
  
  SMCCriterion* smc_criterion = new CESSSMCCriterion(annealing_desired_cess_in);
  //this->model_and_algorithm.smc_termination = NULL;
  this->sequencer = Sequencer(this->the_worker,
                              schedule_in,
                              variable_in,
                              number_of_bisections_in,
                              smc_criterion);
  //this->sequencer_parameters = &this->sequencer.schedule_parameters;
  
  this->mcmc = mcmc_in;
  this->mcmc->set_index(new VectorSingleIndex(indices));
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
:SMC(rng_in, seed_in, data_in, algorithm_parameters, number_of_particles_in, std::max<size_t>(2,lag_in), lag_proposed_in, mcmc_in->get_proposals(), resampling_desired_ess_in, proposal_is_evaluated_in, smcfixed_flag_in, sequencer_limit_is_fixed_in, transform_proposed_particles, results_name_in)
{
  //Parameters candidate_parameters = proposal_in->independent_simulate(*this->rng);
  
  mcmc_in->set_proposal_parameters(&this->algorithm_parameters);
  proposal_in->set_proposal_parameters(&this->algorithm_parameters);
  
  std::vector<size_t> indices;
  indices.reserve(likelihood_estimators_in.size());
  for (size_t i=0; i<likelihood_estimators_in.size(); ++i)
    indices.push_back(i);
  this->index = new VectorSingleIndex(indices);
  
  this->factors = new VectorFactors(likelihood_estimators_in);
  
  /*
  for (auto i=likelihood_estimators.begin();
       i!=likelihood_estimators.end();
       ++i)
  {
    (*i)->setup(candidate_parameters);
  }
  */
  
  // Need to construct LikelihoodEstimator to read in to this constructor.
  this->particle_simulator = new ParameterParticleSimulator(proposal_in,
                                                            likelihood_estimators_in);
  
  if (parallel_in==true)
  {
    //this->the_worker = new RcppParallelSMCWorker(this, grain_size_in);
    Rcpp::stop("Parallel worker not set up.");
  }
  else
  {
    this->the_worker = new SequentialSMCWorker(this);
  }
  
  SMCCriterion* smc_criterion = new CESSSMCCriterion(annealing_desired_cess_in);
  //this->model_and_algorithm.smc_termination = NULL;
  this->sequencer = Sequencer(this->the_worker,
                              schedule_in,
                              sequence_variable_in,
                              number_of_bisections_in,
                              smc_criterion,
                              termination_in);
  //this->sequencer_parameters = &this->sequencer.schedule_parameters;
  
  this->mcmc = mcmc_in;
  if (this->mcmc!=NULL)
    this->mcmc->set_index_if_null(new VectorSingleIndex(indices));
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
:SMC(rng_in, seed_in, data_in, algorithm_parameters, number_of_particles_in, std::max<size_t>(2,lag_in), lag_proposed_in, mcmc_in->get_proposals(), adaptive_resampling_in, proposal_is_evaluated_in, smcfixed_flag_in, sequencer_limit_is_fixed_in, transform_proposed_particles, results_name_in)
{
  //Parameters candidate_parameters = proposal_in->independent_simulate(*this->rng);
  
  mcmc_in->set_proposal_parameters(&this->algorithm_parameters);
  proposal_in->set_proposal_parameters(&this->algorithm_parameters);
  
  //std::vector<size_t> indices;
  //indices.reserve(likelihood_estimators_in.size());
  //for (size_t i=0; i<likelihood_estimators_in.size(); ++i)
  //  indices.push_back(i);
  //this->index = new VectorSingleIndex(indices);
  this->index = without_cancelled_index;
  
  this->factors = new VectorFactors(likelihood_estimators_in);
  
  /*
   for (auto i=likelihood_estimators.begin();
   i!=likelihood_estimators.end();
   ++i)
   {
   (*i)->setup(candidate_parameters);
   }
   */
  
  // Need to construct LikelihoodEstimator to read in to this constructor.
  this->particle_simulator = new ParameterParticleSimulator(proposal_in,
                                                            likelihood_estimators_in);
  
  //Rcout << "Test" << std::endl;
  
  if (parallel_in==true)
  {
    //this->the_worker = new RcppParallelSMCWorker(this,grain_size_in);
    Rcpp::stop("Parallel worker not set up.");
  }
  else
  {
    this->the_worker = new SequentialSMCWorker(this);
  }
  
  SMCCriterion* smc_criterion = adaptive_target_in;
  //this->model_and_algorithm.smc_termination = NULL;
  this->sequencer = Sequencer(this->the_worker,
                              schedules_in,
                              sequence_variables_in,
                              number_of_bisections_in,
                              smc_criterion,
                              termination_in);
  //this->sequencer_parameters = &this->sequencer.schedule_parameters;
  
  this->mcmc = mcmc_in;
  if (this->mcmc!=NULL)
    this->mcmc->set_index_if_null(full_index);
                                  
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
  if (this->mcmc!=NULL)
    delete this->mcmc;
  
  if (this->index!=NULL)
    delete this->index;
}

void SMCMCMCMove::operator=(const SMCMCMCMove &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  if (this->mcmc!=NULL)
    delete this->mcmc;
  
  if (this->index!=NULL)
    delete this->index;

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
    this->mcmc = another.mcmc->mcmc_duplicate();
  else
    this->mcmc = NULL;
  
  if (another.index!=NULL)
    this->index = another.index->duplicate();
  else
    this->index = NULL;
  
  this->mcmc_at_last_step = another.mcmc_at_last_step;
}

SMCOutput* SMCMCMCMove::specific_run()
{
  //if (this->index==NULL)
  //{
  //  Rcout << "Full index is NULL." << std::endl;
  //}
  
  SMCOutput* simulation = this->initialise_smc();
  this->simulate_smc(simulation);
  this->evaluate_smc(simulation);
  simulation->normalise_and_resample_weights();
  return simulation;
}

/*
SMCOutput* SMCMCMCMove::specific_run(const std::string &directory_name)
{
  SMCOutput* simulation = this->initialise_smc();
  this->simulate_smc(simulation);
  this->evaluate_smc(simulation);
  simulation->normalise_and_resample_weights();
  return simulation;
}
*/

SMCOutput* SMCMCMCMove::specific_initialise_smc()
{
  SMCOutput* output = new SMCOutput(this, this->lag, this->lag_proposed, this->results_name);
  return output;
}

void SMCMCMCMove::simulate_smc(SMCOutput* current_state)
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
    //current_state->resample();
    
    //std::ofstream test_file_stream;
    //test_file_stream.open("/Users/richard/Dropbox/code/ilike/experiments/test2.txt",std::ios::out | std::ios::app);
    
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Particles* current_particles = &current_state->back();
    Particles* next_particles = current_state->add_particles(current_particles);
    
    //test_file_stream << "p1" << std::endl;
    
    this->the_worker->move(next_particles,
                           current_particles);
    
    //test_file_stream << "p2" << std::endl;
    
    //this->evaluate_smcfixed_part_smc(current_state);
    current_state->increment_smc_iteration();
    
    //test_file_stream << "p3" << std::endl;
    
    //test_file_stream.close();
    
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
  this->the_worker->smcfixed_weight(this->index,
                                    current_state->back());
  //current_state->initialise_next_step();
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
    
    //std::ofstream test_file_stream;
    //test_file_stream.open("/Users/richard/Dropbox/code/ilike/experiments/test.txt",std::ios::out | std::ios::app);
    
    //test_file_stream << "1" << std::endl;
    
    //Rcout << "p1" << std::endl;
    
    this->sequencer.find_desired_criterion(current_state);
    
    //test_file_stream << "2" << std::endl;
    
    //Rcout << "p2" << std::endl;
    
    // (involves evaluating adaptive weights, using Sequencer)
    this->sequencer.find_next_target_bisection(current_state,this->index);
    
    //test_file_stream << "3" << std::endl;
    
    //Rcout << "p3" << std::endl;

    current_state->log_likelihood = current_state->log_likelihood + current_state->calculate_latest_log_normalising_constant_ratio();
    current_state->llhds.push_back(current_state->log_likelihood);
    
    //test_file_stream << "4" << std::endl;
    
    //Rcout << "p4" << std::endl;
    
    //if (this->sequencer_parameters!=NULL)
    //  current_state->back().schedule_parameters = *this->sequencer_parameters;
    current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();
    
    //test_file_stream << "4" << std::endl;
    
    //Rcout << "p4" << std::endl;
    
    //Rcout << current_state->back().schedule_parameters << std::endl;
    
    //if (current_state)
    
    //this->the_worker->smcadaptive_given_smcfixed_weight(conditioned_on_parameters);
    //current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    // check termination, using sequencer
    if (this->sequencer.check_termination())
    {
      //terminate = TRUE;
      
      if (this->mcmc_at_last_step)
      {
        this->simulate_smc(current_state);
        current_state->decrement_smc_iteration();
      }
      break;
    }
    
    //test_file_stream << "5" << std::endl;
    
    this->simulate_smc(current_state);
    
    //test_file_stream << "6" << std::endl;
    
    //test_file_stream.close();
    
    current_state->back().set_previous_target_evaluated_to_target_evaluated();
  }
}

MoveOutput* SMCMCMCMove::move(RandomNumberGenerator &rng,
                              const Particle &particle) const
{
  return this->mcmc->run(rng,
                         particle);
}

//void SMCMCMCMove::weight_for_adapting_sequence(Particles &current_particles)
//{
//  this->the_worker->smcadaptive_given_smcfixed_weight(this->index,
//                                                      current_particles);
//}

SMCOutput* SMCMCMCMove::specific_run(const Parameters &conditioned_on_parameters)
{
  SMCOutput* simulation = this->initialise_smc(conditioned_on_parameters);
  this->simulate_smc(simulation, conditioned_on_parameters);
  this->evaluate_smc(simulation, conditioned_on_parameters);
  simulation->normalise_and_resample_weights();
  return simulation;
}

/*
SMCOutput* SMCMCMCMove::specific_run(const std::string &directory_name,
                                     const Parameters &conditioned_on_parameters)
{
  SMCOutput* simulation = this->initialise_smc(conditioned_on_parameters);
  this->simulate_smc(simulation, conditioned_on_parameters);
  this->evaluate_smc(simulation, conditioned_on_parameters);
  simulation->normalise_and_resample_weights();
  return simulation;
}
*/

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
    // update MCMC proposals
    this->mcmc->smc_adapt(current_state);
    
    current_state->normalise_and_resample_weights();
    //current_state->resample();
    
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Particles* current_particles = &current_state->back();
    Particles* next_particles = current_state->add_particles(current_particles);
    this->the_worker->move(next_particles,
                           current_particles);
    
    //this->evaluate_smcfixed_part_smc(current_state,
    //                                 conditioned_on_parameters);
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
  //current_state->initialise_next_step();
}

void SMCMCMCMove::evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* current_state,
                                                               const Parameters &conditioned_on_parameters)
{
  // set sequencer to have values from conditioned_on_parameters
  if (!this->sequencer_limit_is_fixed)
    this->sequencer.set_next_with_parameter(conditioned_on_parameters);
  else
    this->sequencer.reset();
  
  // iterate until stop.
  bool terminate = FALSE;
  while (terminate==FALSE)
  {
    this->sequencer.find_desired_criterion(current_state);
    
    // (involves evaluating adaptive weights, using Sequencer)
    this->sequencer.find_next_target_bisection(current_state,
                                               this->index);
    current_state->log_likelihood = current_state->log_likelihood + current_state->calculate_latest_log_normalising_constant_ratio();
    current_state->llhds.push_back(current_state->log_likelihood);
    
    current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();
    
    //this->the_worker->smcadaptive_given_smcfixed_weight(conditioned_on_parameters);
    //current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    // check termination, using sequencer
    if (this->sequencer.check_termination())
    {
      //terminate = TRUE;
      if (this->mcmc_at_last_step)
      {
        this->simulate_smc(current_state,conditioned_on_parameters);
        current_state->decrement_smc_iteration();
      }
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
    //current_state->resample();
    
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Particles* current_particles = &current_state->back();
    Particles* next_particles = current_state->add_particles(current_particles);
    this->the_worker->subsample_move(next_particles,
                                     current_particles);
    
    //this->subsample_evaluate_smcfixed_part_smc(current_state);
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
    //current_state->resample();
    
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Particles* current_particles = &current_state->back();
    Particles* next_particles = current_state->add_particles(current_particles);
    this->the_worker->subsample_move(next_particles,
                                     current_particles);
    
    //this->subsample_evaluate_smcfixed_part_smc(current_state,
    //                                           conditioned_on_parameters);
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
  //current_state->initialise_next_step();
}

void SMCMCMCMove::subsample_evaluate_smcfixed_part_smc(SMCOutput* current_state,
                                                       const Parameters &conditioned_on_parameters)
{
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
  //current_state->initialise_next_step();
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
    this->sequencer.find_desired_criterion(current_state);
    
    // (involves evaluating adaptive weights, using Sequencer)
    this->sequencer.subsample_find_next_target_bisection(current_state,
                                                         this->index);
    
    current_state->subsample_log_likelihood = current_state->subsample_log_likelihood + current_state->calculate_latest_log_normalising_constant_ratio();
    
    //if (this->sequencer_parameters!=NULL)
    //  current_state->back().schedule_parameters = *this->sequencer_parameters;
    current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();
    
    //this->the_worker->smcadaptive_given_smcfixed_weight(conditioned_on_parameters);
    //current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    // check termination, using sequencer
    if (this->sequencer.check_termination())
    {
      //terminate = TRUE;
      if (this->mcmc_at_last_step)
      {
        this->simulate_smc(current_state);
        current_state->decrement_smc_iteration();
      }
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
    this->sequencer.find_desired_criterion(current_state);
    
    // (involves evaluating adaptive weights, using Sequencer)
    this->sequencer.subsample_find_next_target_bisection(current_state,
                                                         this->index);
    
    current_state->subsample_log_likelihood = current_state->subsample_log_likelihood + current_state->calculate_latest_log_normalising_constant_ratio();
    
    //if (this->sequencer_parameters!=NULL)
    //  current_state->back().schedule_parameters = *this->sequencer_parameters;
    current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();
    
    //this->the_worker->smcadaptive_given_smcfixed_weight(conditioned_on_parameters);
    //current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    // check termination, using sequencer
    if (this->sequencer.check_termination())
    {
      //terminate = TRUE;
      if (this->mcmc_at_last_step)
      {
        this->simulate_smc(current_state,conditioned_on_parameters);
        current_state->decrement_smc_iteration();
      }
      break;
    }
    
    this->subsample_simulate_smc(current_state,
                                 conditioned_on_parameters);
    
    current_state->back().subsample_set_previous_target_evaluated_to_target_evaluated();
  }
}

/*
MoveOutput* SMCMCMCMove::move(RandomNumberGenerator &rng,
                              Particle &particle,
                              const Parameters &conditioned_on_parameters)
{
  return this->mcmc->run(rng,
                         particle,
                         conditioned_on_parameters);
}
*/

void SMCMCMCMove::weight_for_adapting_sequence(const Index* index,
                                               Particles &current_particles)
{
  this->the_worker->smcadaptive_given_smcfixed_weight(this->index,
                                                      current_particles);
}

/*
void SMCMCMCMove::weight_for_adapting_sequence(const Index* index,
                                               Particles &current_particles,
                                               const Parameters &conditioned_on_parameters)
{
  this->the_worker->smcadaptive_given_smcfixed_weight(this->index,
                                                      current_particles,
                                                      conditioned_on_parameters);
}
*/

MoveOutput* SMCMCMCMove::subsample_move(RandomNumberGenerator &rng,
                                        const Particle &particle) const
{
  return this->mcmc->subsample_run(rng,
                                   particle);
}

/*
MoveOutput* SMCMCMCMove::subsample_move(RandomNumberGenerator &rng,
                                        Particle &particle,
                                        const Parameters &conditioned_on_parameters)
{
  return this->mcmc->subsample_run(rng,
                                   particle,
                                   conditioned_on_parameters);
}
*/

void SMCMCMCMove::subsample_weight_for_adapting_sequence(const Index* index,
                                                         Particles &current_particles)
{
  this->the_worker->subsample_smcadaptive_given_smcfixed_weight(this->index,
                                                                current_particles);
}

/*
void SMCMCMCMove::subsample_weight_for_adapting_sequence(const Index* index,
                                                         Particles &current_particles,
                                                         const Parameters &conditioned_on_parameters)
{
  this->the_worker->subsample_smcadaptive_given_smcfixed_weight(this->index,
                                                                current_particles,
                                                                conditioned_on_parameters);
}
*/

// void SMCMCMCMove::smc_step(void)
// {
//
// }
//
// void SMCMCMCMove::weight_update(void)
// {
//
// }

//void SMCMCMCMove::smc_update(SMCOutput* current_state)
//{
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

//}
