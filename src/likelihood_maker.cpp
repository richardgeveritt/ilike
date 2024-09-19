#include "ilike_header.h"
#include "exact_likelihood_estimator.h"
#include "abc_likelihood_estimator.h"
#include "lp_uniform_abc_kernel_factor.h"
#include "gaussian_abc_kernel_factor.h"
#include "custom_distribution_proposal_kernel.h"
#include "custom_independent_proposal_kernel.h"
#include "transformed_independent_proposal_kernel.h"
#include "importance_sampler.h"
#include "smc_mcmc_move.h"
#include "always_smc_termination.h"
#include "density_likelihood_estimator.h"
#include "gaussian_density_estimator.h"
#include "annealed_likelihood_estimator.h"
#include "utils.h"
#include "stable_smc_termination.h"
#include "vector_index.h"
#include "ensemble_kalman_inversion.h"
#include "metropolis_mcmc.h"
#include "ensemble_kalman_filter.h"

namespace ilike
{

ABCLikelihoodEstimator* make_lp_uniform_abc_kernel(RandomNumberGenerator* rng_in,
                                                   size_t* seed_in,
                                                   Data* data_in,
                                                   double p_in,
                                                   const std::vector<std::string> &data_variables_in,
                                                   const std::string &scale_variable_in,
                                                   const std::string &epsilon_variable_in,
                                                   bool fixed_epsilon_in)
{
  return new ABCLikelihoodEstimator(rng_in,
                                    seed_in,
                                    data_in,
                                    new LpUniformABCKernelFactor(p_in,
                                                                 data_variables_in,
                                                                 epsilon_variable_in,
                                                                 scale_variable_in,
                                                                 data_in),
                                    fixed_epsilon_in);
}

ABCLikelihoodEstimator* make_gaussian_abc_kernel(RandomNumberGenerator* rng_in,
                                                 size_t* seed_in,
                                                 Data* data_in,
                                                 const std::vector<std::string> &data_variables_in,
                                                 const std::string &scale_variable_in,
                                                 const std::string &epsilon_variable_in,
                                                 bool fixed_epsilon_in)
{
  return new ABCLikelihoodEstimator(rng_in,
                                    seed_in,
                                    data_in,
                                    new GaussianABCKernelFactor(data_variables_in,
                                                                epsilon_variable_in,
                                                                scale_variable_in,
                                                                data_in),
                                    fixed_epsilon_in);
}

ImportanceSampler* make_fixed_epsilon_lp_uniform_abc_likelihood(RandomNumberGenerator* rng_in,
                                                                size_t* seed_in,
                                                                Data* data_in,
                                                                double p_in,
                                                                const std::vector<std::string> &data_variables_in,
                                                                const std::string &scale_variable_in,
                                                                const std::string &epsilon_variable_in,
                                                                double epsilon_in,
                                                                SimulateModelPtr simulate_model_in,
                                                                size_t number_of_abc_simulations_in,
                                                                bool store_output_in,
                                                                bool parallel_in,
                                                                size_t grain_size_in)
{
  std::vector<LikelihoodEstimator*> abc_kernel;
  abc_kernel.push_back(make_lp_uniform_abc_kernel(rng_in,
                                                  seed_in,
                                                  data_in,
                                                  p_in,
                                                  data_variables_in,
                                                  scale_variable_in,
                                                  epsilon_variable_in,
                                                  true));
  
  IndependentProposalKernel* model_simulator = new CustomIndependentProposalKernel(simulate_model_in);
  
  return new ImportanceSampler(rng_in,
                               seed_in,
                               data_in,
                               Parameters(),
                               number_of_abc_simulations_in,
                               epsilon_variable_in,
                               epsilon_in,
                               abc_kernel,
                               model_simulator,
                               false,
                               store_output_in,
                               true,
                               true,
                               false,
                               parallel_in,
                               grain_size_in,
                               "");
}

ImportanceSampler* make_fixed_epsilon_lp_uniform_abc_likelihood(RandomNumberGenerator* rng_in,
                                                                size_t* seed_in,
                                                                Data* summary_data_in,
                                                                double p_in,
                                                                const std::vector<std::string> &summary_data_variables_in,
                                                                const std::string &scale_variable_in,
                                                                const std::string &epsilon_variable_in,
                                                                double epsilon_in,
                                                                SimulateModelPtr simulate_model_in,
                                                                std::shared_ptr<Transform> summary_statistics_in,
                                                                size_t number_of_abc_simulations_in,
                                                                bool store_output_in,
                                                                bool parallel_in,
                                                                size_t grain_size_in)
{
  std::vector<LikelihoodEstimator*> abc_kernel;
  abc_kernel.push_back(make_lp_uniform_abc_kernel(rng_in,
                                                  seed_in,
                                                  summary_data_in,
                                                  p_in,
                                                  summary_data_variables_in,
                                                  scale_variable_in,
                                                  epsilon_variable_in,
                                                  true));
  
  IndependentProposalKernel* model_simulator = new CustomIndependentProposalKernel(simulate_model_in);
  IndependentProposalKernel* summary_model_simulator = new TransformedIndependentProposalKernel(model_simulator,
                                                                                                summary_statistics_in,
                                                                                                true);
  
  return new ImportanceSampler(rng_in,
                               seed_in,
                               summary_data_in,
                               Parameters(),
                               number_of_abc_simulations_in,
                               epsilon_variable_in,
                               epsilon_in,
                               abc_kernel,
                               summary_model_simulator,
                               false,
                               store_output_in,
                               true,
                               true,
                               false,
                               parallel_in,
                               grain_size_in,
                               "");
}

ImportanceSampler* make_varying_epsilon_lp_uniform_abc_likelihood(RandomNumberGenerator* rng_in,
                                                                  size_t* seed_in,
                                                                  Data* data_in,
                                                                  double p_in,
                                                                  const std::vector<std::string> &data_variables_in,
                                                                  const std::string &scale_variable_in,
                                                                  const std::string &epsilon_variable_in,
                                                                  SimulateModelPtr simulate_model_in,
                                                                  size_t number_of_abc_simulations_in,
                                                                  bool store_output_in,
                                                                  bool parallel_in,
                                                                  size_t grain_size_in)
{
  std::vector<LikelihoodEstimator*> abc_kernel;
  abc_kernel.push_back(make_lp_uniform_abc_kernel(rng_in,
                                                  seed_in,
                                                  data_in,
                                                  p_in,
                                                  data_variables_in,
                                                  scale_variable_in,
                                                  epsilon_variable_in,
                                                  false));
  
  IndependentProposalKernel* model_simulator = new CustomIndependentProposalKernel(simulate_model_in);
  
  return new ImportanceSampler(rng_in,
                               seed_in,
                               data_in,
                               Parameters(),
                               number_of_abc_simulations_in,
                               epsilon_variable_in,
                               abc_kernel,
                               model_simulator,
                               false,
                               store_output_in,
                               false,
                               true,
                               false,
                               parallel_in,
                               grain_size_in,
                               "");
}

ImportanceSampler* make_varying_epsilon_lp_uniform_abc_likelihood(RandomNumberGenerator* rng_in,
                                                                  size_t* seed_in,
                                                                  Data* summary_data_in,
                                                                  double p_in,
                                                                  const std::vector<std::string> &summary_data_variables_in,
                                                                  const std::string &scale_variable_in,
                                                                  const std::string &epsilon_variable_in,
                                                                  SimulateModelPtr simulate_model_in,
                                                                  std::shared_ptr<Transform> summary_statistics_in,
                                                                  size_t number_of_abc_simulations_in,
                                                                  bool store_output_in,
                                                                  bool parallel_in,
                                                                  size_t grain_size_in)
{
  std::vector<LikelihoodEstimator*> abc_kernel;
  abc_kernel.push_back(make_lp_uniform_abc_kernel(rng_in,
                                                  seed_in,
                                                  summary_data_in,
                                                  p_in,
                                                  summary_data_variables_in,
                                                  scale_variable_in,
                                                  epsilon_variable_in,
                                                  false));
  
  IndependentProposalKernel* model_simulator = new CustomIndependentProposalKernel(simulate_model_in);
  IndependentProposalKernel* summary_model_simulator = new TransformedIndependentProposalKernel(model_simulator,
                                                                                                summary_statistics_in,
                                                                                                false);
  
  return new ImportanceSampler(rng_in,
                               seed_in,
                               summary_data_in,
                               Parameters(),
                               number_of_abc_simulations_in,
                               epsilon_variable_in,
                               abc_kernel,
                               summary_model_simulator,
                               false,
                               store_output_in,
                               false,
                               false,
                               false,
                               parallel_in,
                               grain_size_in,
                               "");
}



ImportanceSampler* make_fixed_epsilon_gaussian_abc_likelihood(RandomNumberGenerator* rng_in,
                                                              size_t* seed_in,
                                                              Data* data_in,
                                                              const std::vector<std::string> &data_variables_in,
                                                              const std::string &scale_variable_in,
                                                              const std::string &epsilon_variable_in,
                                                              double epsilon_in,
                                                              SimulateModelPtr simulate_model_in,
                                                              size_t number_of_abc_simulations_in,
                                                              bool store_output_in,
                                                              bool parallel_in,
                                                              size_t grain_size_in)
{
  std::vector<LikelihoodEstimator*> abc_kernel;
  abc_kernel.push_back(make_gaussian_abc_kernel(rng_in,
                                                seed_in,
                                                data_in,
                                                data_variables_in,
                                                scale_variable_in,
                                                epsilon_variable_in,
                                                true));
  
  IndependentProposalKernel* model_simulator = new CustomIndependentProposalKernel(simulate_model_in);
  
  return new ImportanceSampler(rng_in,
                               seed_in,
                               data_in,
                               Parameters(),
                               number_of_abc_simulations_in,
                               epsilon_variable_in,
                               epsilon_in,
                               abc_kernel,
                               model_simulator,
                               false,
                               store_output_in,
                               true,
                               true,
                               false,
                               parallel_in,
                               grain_size_in,
                               "");
}

EnsembleKalmanInversion* make_fixed_epsilon_enki_abc_likelihood(RandomNumberGenerator* rng_in,
                                                                size_t* seed_in,
                                                                Data* data_in,
                                                                size_t lag_in,
                                                                EnsembleShifter* shifter_in,
                                                                //double annealing_desired_cess_in,
                                                                //size_t number_of_bisections_in,
                                                                const std::string &scale_variable_in,
                                                                const std::string &sequence_variable_in,
                                                                //const std::vector<double> &schedule_in,
                                                                double min_epsilon,
                                                                size_t number_of_targets_in,
                                                                const arma::colvec &scale_in,
                                                                const std::vector<std::string> &data_variables_in,
                                                                SimulateModelPtr simulate_model_in,
                                                                size_t number_of_abc_simulations_in,
                                                                double significance_level_in,
                                                                size_t estimator_type_in,
                                                                bool parallel_in,
                                                                size_t grain_size_in)
{
  IndependentProposalKernel* model_simulator = new CustomIndependentProposalKernel(simulate_model_in);
  
  /*
   if (schedule_in.size()<2)
   {
   stop("make_fixed_epsilon_enki_abc_likelihood - schedule needs to be of length at least 2.");
   }
   double min_epsilon = schedule_in[0];
   for (size_t i=1; i<schedule_in.size(); ++i)
   {
   if (schedule_in[i]<schedule_in[i-1])
   {
   min_epsilon = schedule_in[i];
   }
   else
   {
   stop("make_fixed_epsilon_enki_abc_likelihood - schedule needs to be strictly decreasing.");
   }
   }
   
   std::vector<double> temperature_schedule;
   temperature_schedule.reserve(schedule_in.size());
   for (size_t i=0; i<schedule_in.size(); ++i)
   {
   if (schedule_in[i]==arma::datum::inf)
   {
   temperature_schedule.push_back(0.0);
   }
   else
   {
   temperature_schedule.push_back(pow(min_epsilon/schedule_in[i],2.0));
   }
   }
   */
  
  if (min_epsilon<=0.0)
  {
    stop("make_fixed_epsilon_enki_abc_likelihood - tolerance must be positive.");
  }
  
  return new EnsembleKalmanInversion(rng_in,
                                     seed_in,
                                     data_in,
                                     number_of_abc_simulations_in,
                                     lag_in,
                                     shifter_in,
                                     //annealing_desired_cess_in,
                                     number_of_targets_in,
                                     //number_of_bisections_in,
                                     sequence_variable_in,
                                     //temperature_schedule,
                                     model_simulator,
                                     data_variables_in,
                                     min_epsilon,
                                     scale_variable_in,
                                     scale_in,
                                     NULL,
                                     NULL,
                                     significance_level_in,
                                     estimator_type_in,
                                     parallel_in,
                                     grain_size_in,
                                     "");
}

EnsembleKalmanFilter* make_enkf_likelihood(RandomNumberGenerator* rng_in,
                                           size_t* seed_in,
                                           Data* data_in,
                                           size_t lag_in,
                                           const std::string &index_name_in,
                                           const std::string &time_name_in,
                                           const std::string &time_diff_name_in,
                                           size_t first_index_in,
                                           size_t last_index_in,
                                           size_t predictions_per_update_in,
                                           double update_time_step_in,
                                           double initial_time_in,
                                           size_t number_of_ensemble_members_in,
                                           EnsembleShifter* shifter_in,
                                           IndependentProposalKernel* prior_in,
                                           ProposalKernel* transition_model_in,
                                           const std::vector<MeasurementCovarianceEstimator*> &measurement_covariance_estimators_in,
                                           bool parallel_in,
                                           size_t grain_size_in)
{
  return new EnsembleKalmanFilter(rng_in,
                                  seed_in,
                                  data_in,
                                  lag_in,
                                  index_name_in,
                                  time_name_in,
                                  time_diff_name_in,
                                  first_index_in,
                                  last_index_in,
                                  predictions_per_update_in,
                                  update_time_step_in,
                                  initial_time_in,
                                  number_of_ensemble_members_in,
                                  shifter_in,
                                  NULL,
                                  prior_in,
                                  transition_model_in,
                                  measurement_covariance_estimators_in,
                                  TRUE,
                                  TRUE,
                                  parallel_in,
                                  grain_size_in,
                                  "");
}

/*
 SMCMCMCMove* make_fixed_epsilon_smc_gaussian_abc_likelihood(RandomNumberGenerator* rng_in,
 size_t* seed_in,
 Data* data_in,
 size_t lag_in,
 MCMC* mcmc_in,
 SMCCriterion* adaptive_resampling_in,
 SMCCriterion* adaptive_target_in,
 size_t number_of_bisections_in,
 const std::string &scale_variable_in,
 const std::string &sequence_variable_in,
 const std::vector<double> &schedule_in,
 const std::vector<std::string> &data_variables_in,
 SimulateModelPtr simulate_model_in,
 size_t number_of_abc_simulations_in,
 bool parallel_in,
 size_t grain_size_in)
 {
 std::vector<LikelihoodEstimator*> abc_kernel;
 abc_kernel.push_back(make_gaussian_abc_kernel(rng_in,
 seed_in,
 data_in,
 data_variables_in,
 scale_variable_in,
 epsilon_variable_in,
 true));
 
 std::vector<std::string> sequence_variables_in;
 sequence_variables_in.push_back(sequence_variable_in);
 
 std::vector<std::vector<double>> schedules_in
 
 return new SMCMCMCMove(rng_in,
 seed_in,
 data_in,
 Parameters(),
 number_of_abc_simulations_in,
 lag_in,
 lag_in,
 mcmc_in,
 adaptive_resampling_in,
 adaptive_target_in,
 number_of_bisections_in,
 NULL,
 sequence_variables_in,
 schedules_in,
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
 const std::string &results_name_in);
 }
 */

ImportanceSampler* make_fixed_epsilon_gaussian_abc_likelihood(RandomNumberGenerator* rng_in,
                                                              size_t* seed_in,
                                                              Data* summary_data_in,
                                                              const std::vector<std::string> &summary_data_variables_in,
                                                              const std::string &scale_variable_in,
                                                              const std::string &epsilon_variable_in,
                                                              double epsilon_in,
                                                              SimulateModelPtr simulate_model_in,
                                                              std::shared_ptr<Transform> summary_statistics_in,
                                                              size_t number_of_abc_simulations_in,
                                                              bool store_output_in,
                                                              bool parallel_in,
                                                              size_t grain_size_in)
{
  std::vector<LikelihoodEstimator*> abc_kernel;
  abc_kernel.push_back(make_gaussian_abc_kernel(rng_in,
                                                seed_in,
                                                summary_data_in,
                                                summary_data_variables_in,
                                                scale_variable_in,
                                                epsilon_variable_in,
                                                true));
  
  IndependentProposalKernel* model_simulator = new CustomIndependentProposalKernel(simulate_model_in);
  IndependentProposalKernel* summary_model_simulator = new TransformedIndependentProposalKernel(model_simulator,
                                                                                                summary_statistics_in,
                                                                                                true);
  
  return new ImportanceSampler(rng_in,
                               seed_in,
                               summary_data_in,
                               Parameters(),
                               number_of_abc_simulations_in,
                               epsilon_variable_in,
                               epsilon_in,
                               abc_kernel,
                               summary_model_simulator,
                               false,
                               store_output_in,
                               true,
                               true,
                               false,
                               parallel_in,
                               grain_size_in,
                               "");
}

EnsembleKalmanInversion* make_fixed_epsilon_enki_abc_likelihood(RandomNumberGenerator* rng_in,
                                                                size_t* seed_in,
                                                                Data* summary_data_in,
                                                                size_t lag_in,
                                                                EnsembleShifter* shifter_in,
                                                                //double annealing_desired_cess_in,
                                                                //size_t number_of_bisections_in,
                                                                const std::string &scale_variable_in,
                                                                const std::string &sequence_variable_in,
                                                                //const std::vector<double> &schedule_in,
                                                                double min_epsilon,
                                                                size_t number_of_targets_in,
                                                                const arma::colvec &scale_in,
                                                                const std::vector<std::string> &summary_data_variables_in,
                                                                SimulateModelPtr simulate_model_in,
                                                                std::shared_ptr<Transform> summary_statistics_in,
                                                                bool enki_on_summary,
                                                                size_t number_of_abc_simulations_in,
                                                                double significance_level_in,
                                                                size_t estimator_type_in,
                                                                bool parallel_in,
                                                                size_t grain_size_in)
{
  IndependentProposalKernel* model_simulator = new CustomIndependentProposalKernel(simulate_model_in);
  if (enki_on_summary)
  {
    model_simulator = new TransformedIndependentProposalKernel(model_simulator,
                                                               summary_statistics_in);
  }
  
  /*
   if (schedule_in.size()<2)
   {
   stop("make_fixed_epsilon_enki_abc_likelihood - schedule needs to be of length at least 2.");
   }
   double min_epsilon = schedule_in[0];
   for (size_t i=1; i<schedule_in.size(); ++i)
   {
   if (schedule_in[i]<schedule_in[i-1])
   {
   min_epsilon = schedule_in[i];
   }
   else
   {
   stop("make_fixed_epsilon_enki_abc_likelihood - schedule needs to be strictly decreasing.");
   }
   }
   
   std::vector<double> temperature_schedule;
   temperature_schedule.reserve(schedule_in.size());
   for (size_t i=0; i<schedule_in.size(); ++i)
   {
   if (schedule_in[i]==arma::datum::inf)
   {
   temperature_schedule.push_back(0.0);
   }
   else
   {
   temperature_schedule.push_back(pow(min_epsilon/schedule_in[i],2.0));
   }
   }
   */
  
  if (min_epsilon<=0.0)
  {
    stop("make_fixed_epsilon_enki_abc_likelihood - tolerance must be positive.");
  }
  
  if (enki_on_summary)
  {
    return new EnsembleKalmanInversion(rng_in,
                                       seed_in,
                                       summary_data_in,
                                       number_of_abc_simulations_in,
                                       lag_in,
                                       shifter_in,
                                       //annealing_desired_cess_in,
                                       number_of_targets_in,
                                       //number_of_bisections_in,
                                       sequence_variable_in,
                                       //temperature_schedule,
                                       model_simulator,
                                       summary_data_variables_in,
                                       min_epsilon,
                                       scale_variable_in,
                                       scale_in,
                                       NULL,
                                       NULL,
                                       significance_level_in,
                                       estimator_type_in,
                                       parallel_in,
                                       grain_size_in,
                                       "");
  }
  else
  {
    return new EnsembleKalmanInversion(rng_in,
                                       seed_in,
                                       summary_data_in,
                                       number_of_abc_simulations_in,
                                       lag_in,
                                       shifter_in,
                                       //annealing_desired_cess_in,
                                       number_of_targets_in,
                                       //number_of_bisections_in,
                                       sequence_variable_in,
                                       //temperature_schedule,
                                       model_simulator,
                                       summary_data_variables_in,
                                       min_epsilon,
                                       scale_variable_in,
                                       scale_in,
                                       summary_statistics_in,
                                       NULL,
                                       significance_level_in,
                                       estimator_type_in,
                                       parallel_in,
                                       grain_size_in,
                                       "");
  }
}

ImportanceSampler* make_varying_epsilon_gaussian_abc_likelihood(RandomNumberGenerator* rng_in,
                                                                size_t* seed_in,
                                                                Data* data_in,
                                                                const std::vector<std::string> &data_variables_in,
                                                                const std::string &scale_variable_in,
                                                                const std::string &epsilon_variable_in,
                                                                SimulateModelPtr simulate_model_in,
                                                                size_t number_of_abc_simulations_in,
                                                                bool store_output_in,
                                                                bool parallel_in,
                                                                size_t grain_size_in)
{
  std::vector<LikelihoodEstimator*> abc_kernel;
  abc_kernel.push_back(make_gaussian_abc_kernel(rng_in,
                                                seed_in,
                                                data_in,
                                                data_variables_in,
                                                scale_variable_in,
                                                epsilon_variable_in,
                                                false));
  
  IndependentProposalKernel* model_simulator = new CustomIndependentProposalKernel(simulate_model_in);
  
  return new ImportanceSampler(rng_in,
                               seed_in,
                               data_in,
                               Parameters(),
                               number_of_abc_simulations_in,
                               epsilon_variable_in,
                               abc_kernel,
                               model_simulator,
                               false,
                               store_output_in,
                               false,
                               true,
                               false,
                               parallel_in,
                               grain_size_in,
                               "");
}

ImportanceSampler* make_varying_epsilon_gaussian_abc_likelihood(RandomNumberGenerator* rng_in,
                                                                size_t* seed_in,
                                                                Data* summary_data_in,
                                                                const std::vector<std::string> &summary_data_variables_in,
                                                                const std::string &scale_variable_in,
                                                                const std::string &epsilon_variable_in,
                                                                SimulateModelPtr simulate_model_in,
                                                                std::shared_ptr<Transform> summary_statistics_in,
                                                                size_t number_of_abc_simulations_in,
                                                                bool store_output_in,
                                                                bool parallel_in,
                                                                size_t grain_size_in)
{
  std::vector<LikelihoodEstimator*> abc_kernel;
  abc_kernel.push_back(make_gaussian_abc_kernel(rng_in,
                                                seed_in,
                                                summary_data_in,
                                                summary_data_variables_in,
                                                scale_variable_in,
                                                epsilon_variable_in,
                                                false));
  
  IndependentProposalKernel* model_simulator = new CustomIndependentProposalKernel(simulate_model_in);
  IndependentProposalKernel* summary_model_simulator = new TransformedIndependentProposalKernel(model_simulator,
                                                                                                summary_statistics_in,
                                                                                                false);
  
  return new ImportanceSampler(rng_in,
                               seed_in,
                               summary_data_in,
                               Parameters(),
                               number_of_abc_simulations_in,
                               epsilon_variable_in,
                               abc_kernel,
                               summary_model_simulator,
                               false,
                               store_output_in,
                               false,
                               false,
                               false,
                               parallel_in,
                               grain_size_in,
                               "");
}



DensityLikelihoodEstimator* make_sl_likelihood(RandomNumberGenerator* rng_in,
                                               size_t* seed_in,
                                               Data* data_in,
                                               const std::vector<std::string> &data_variables_in,
                                               SimulateModelPtr simulate_model_in,
                                               size_t number_of_sl_simulations_in,
                                               bool unbiased_in,
                                               bool store_output_in,
                                               bool parallel_in,
                                               size_t grain_size_in)
{
  IndependentProposalKernel* model_simulator = new CustomIndependentProposalKernel(simulate_model_in);
  
  DensityEstimator* density_estimator = new GaussianDensityEstimator(data_variables_in,
                                                                     unbiased_in);
  
  return new DensityLikelihoodEstimator(rng_in,
                                        seed_in,
                                        data_in,
                                        Parameters(),
                                        number_of_sl_simulations_in,
                                        true,
                                        density_estimator,
                                        model_simulator,
                                        false,
                                        store_output_in,
                                        parallel_in,
                                        grain_size_in);
}

DensityLikelihoodEstimator* make_sl_likelihood(RandomNumberGenerator* rng_in,
                                               size_t* seed_in,
                                               Data* summary_data_in,
                                               const std::vector<std::string> &summary_data_variables_in,
                                               SimulateModelPtr simulate_model_in,
                                               std::shared_ptr<Transform> summary_statistics_in,
                                               size_t number_of_sl_simulations_in,
                                               bool unbiased_in,
                                               bool store_output_in,
                                               bool parallel_in,
                                               size_t grain_size_in)
{
  IndependentProposalKernel* model_simulator = new CustomIndependentProposalKernel(simulate_model_in);
  IndependentProposalKernel* summary_model_simulator = new TransformedIndependentProposalKernel(model_simulator,
                                                                                                summary_statistics_in);
  
  DensityEstimator* density_estimator = new GaussianDensityEstimator(summary_data_variables_in,
                                                                     unbiased_in);
  
  return new DensityLikelihoodEstimator(rng_in,
                                        seed_in,
                                        summary_data_in,
                                        Parameters(),
                                        number_of_sl_simulations_in,
                                        true,
                                        density_estimator,
                                        summary_model_simulator,
                                        false,
                                        store_output_in,
                                        parallel_in,
                                        grain_size_in);
}

ImportanceSampler* make_fixed_epsilon_lp_uniform_abc_is(RandomNumberGenerator* rng_in,
                                                        size_t* seed_in,
                                                        Data* data_in,
                                                        size_t number_of_particles_in,
                                                        IndependentProposalKernel* prior_in,
                                                        double p_in,
                                                        const std::vector<std::string> &data_variables_in,
                                                        const std::string &scale_variable_in,
                                                        const std::string &epsilon_variable_in,
                                                        double epsilon_in,
                                                        SimulateModelPtr simulate_model_in,
                                                        size_t number_of_abc_simulations_in,
                                                        bool abc_store_output_in,
                                                        bool abc_parallel_in,
                                                        size_t abc_grain_size_in,
                                                        bool is_store_output_in,
                                                        bool is_parallel_in,
                                                        size_t is_grain_size_in,
                                                        const std::string &results_name_in)
{
  std::vector<LikelihoodEstimator*> abc_likelihood;
  abc_likelihood.push_back(make_fixed_epsilon_lp_uniform_abc_likelihood(rng_in,
                                                                        seed_in,
                                                                        data_in,
                                                                        p_in,
                                                                        data_variables_in,
                                                                        scale_variable_in,
                                                                        epsilon_variable_in,
                                                                        epsilon_in,
                                                                        simulate_model_in,
                                                                        number_of_abc_simulations_in,
                                                                        abc_store_output_in,
                                                                        abc_parallel_in,
                                                                        abc_grain_size_in));
  
  return new ImportanceSampler(rng_in,
                               seed_in,
                               data_in,
                               Parameters(),
                               number_of_particles_in,
                               epsilon_variable_in,
                               epsilon_in,
                               abc_likelihood,
                               prior_in,
                               false,
                               is_store_output_in,
                               true,
                               true,
                               false,
                               is_parallel_in,
                               is_grain_size_in,
                               results_name_in);
}

ImportanceSampler* make_fixed_epsilon_lp_uniform_abc_is(RandomNumberGenerator* rng_in,
                                                        size_t* seed_in,
                                                        Data* summary_data_in,
                                                        size_t number_of_particles_in,
                                                        IndependentProposalKernel* prior_in,
                                                        double p_in,
                                                        const std::vector<std::string> &summary_data_variables_in,
                                                        const std::string &scale_variable_in,
                                                        const std::string &epsilon_variable_in,
                                                        double epsilon_in,
                                                        SimulateModelPtr simulate_model_in,
                                                        std::shared_ptr<Transform> summary_statistics_in,
                                                        size_t number_of_abc_simulations_in,
                                                        bool abc_store_output_in,
                                                        bool abc_parallel_in,
                                                        size_t abc_grain_size_in,
                                                        bool is_store_output_in,
                                                        bool is_parallel_in,
                                                        size_t is_grain_size_in,
                                                        const std::string &results_name_in)
{
  std::vector<LikelihoodEstimator*> abc_likelihood;
  abc_likelihood.push_back(make_fixed_epsilon_lp_uniform_abc_likelihood(rng_in,
                                                                        seed_in,
                                                                        summary_data_in,
                                                                        p_in,
                                                                        summary_data_variables_in,
                                                                        scale_variable_in,
                                                                        epsilon_variable_in,
                                                                        epsilon_in,
                                                                        simulate_model_in,
                                                                        summary_statistics_in,
                                                                        number_of_abc_simulations_in,
                                                                        abc_store_output_in,
                                                                        abc_parallel_in,
                                                                        abc_grain_size_in));
  
  return new ImportanceSampler(rng_in,
                               seed_in,
                               summary_data_in,
                               Parameters(),
                               number_of_particles_in,
                               epsilon_variable_in,
                               epsilon_in,
                               abc_likelihood,
                               prior_in,
                               false,
                               is_store_output_in,
                               true,
                               true,
                               false,
                               is_parallel_in,
                               is_grain_size_in,
                               results_name_in);
}

ImportanceSampler* make_fixed_epsilon_gaussian_abc_is(RandomNumberGenerator* rng_in,
                                                      size_t* seed_in,
                                                      Data* data_in,
                                                      size_t number_of_particles_in,
                                                      IndependentProposalKernel* prior_in,
                                                      const std::vector<std::string> &data_variables_in,
                                                      const std::string &scale_variable_in,
                                                      const std::string &epsilon_variable_in,
                                                      double epsilon_in,
                                                      SimulateModelPtr simulate_model_in,
                                                      size_t number_of_abc_simulations_in,
                                                      bool abc_store_output_in,
                                                      bool abc_parallel_in,
                                                      size_t abc_grain_size_in,
                                                      bool is_store_output_in,
                                                      bool is_parallel_in,
                                                      size_t is_grain_size_in,
                                                      const std::string &results_name_in)
{
  std::vector<LikelihoodEstimator*> abc_likelihood;
  abc_likelihood.push_back(make_fixed_epsilon_gaussian_abc_likelihood(rng_in,
                                                                      seed_in,
                                                                      data_in,
                                                                      data_variables_in,
                                                                      scale_variable_in,
                                                                      epsilon_variable_in,
                                                                      epsilon_in,
                                                                      simulate_model_in,
                                                                      number_of_abc_simulations_in,
                                                                      abc_store_output_in,
                                                                      abc_parallel_in,
                                                                      abc_grain_size_in));
  
  return new ImportanceSampler(rng_in,
                               seed_in,
                               data_in,
                               Parameters(),
                               number_of_particles_in,
                               epsilon_variable_in,
                               epsilon_in,
                               abc_likelihood,
                               prior_in,
                               false,
                               is_store_output_in,
                               true,
                               true,
                               false,
                               is_parallel_in,
                               is_grain_size_in,
                               results_name_in);
}

ImportanceSampler* make_fixed_epsilon_enki_abc_is(RandomNumberGenerator* rng_in,
                                                  size_t* seed_in,
                                                  Data* data_in,
                                                  size_t number_of_particles_in,
                                                  IndependentProposalKernel* prior_in,
                                                  const std::vector<std::string> &data_variables_in,
                                                  const std::string &scale_variable_in,
                                                  const std::string &epsilon_variable_in,
                                                  double epsilon_in,
                                                  size_t number_of_targets_in,
                                                  const arma::colvec &scale_in,
                                                  size_t enki_lag_in,
                                                  EnsembleShifter* shifter_in,
                                                  //double enki_annealing_desired_cess_in,
                                                  //size_t enki_number_of_bisections_in,
                                                  SimulateModelPtr simulate_model_in,
                                                  size_t number_of_abc_simulations_in,
                                                  double significance_level_in,
                                                  size_t estimator_type_in,
                                                  bool abc_parallel_in,
                                                  size_t abc_grain_size_in,
                                                  bool is_store_output_in,
                                                  bool is_parallel_in,
                                                  size_t is_grain_size_in,
                                                  const std::string &results_name_in)
{
  /*
   std::vector<double> enki_schedule;
   enki_schedule.push_back(arma::datum::inf);
   enki_schedule.push_back(epsilon_in);
   */
  
  std::vector<LikelihoodEstimator*> abc_likelihood;
  abc_likelihood.push_back(make_fixed_epsilon_enki_abc_likelihood(rng_in,
                                                                  seed_in,
                                                                  data_in,
                                                                  enki_lag_in,
                                                                  shifter_in,
                                                                  //enki_annealing_desired_cess_in,
                                                                  //enki_number_of_bisections_in,
                                                                  scale_variable_in,
                                                                  epsilon_variable_in,
                                                                  epsilon_in,
                                                                  number_of_targets_in,
                                                                  scale_in,
                                                                  //enki_schedule,
                                                                  data_variables_in,
                                                                  simulate_model_in,
                                                                  number_of_abc_simulations_in,
                                                                  significance_level_in,
                                                                  estimator_type_in,
                                                                  abc_parallel_in,
                                                                  abc_grain_size_in));
  
  return new ImportanceSampler(rng_in,
                               seed_in,
                               data_in,
                               Parameters(),
                               number_of_particles_in,
                               epsilon_variable_in,
                               epsilon_in,
                               abc_likelihood,
                               prior_in,
                               false,
                               is_store_output_in,
                               true,
                               true,
                               false,
                               is_parallel_in,
                               is_grain_size_in,
                               results_name_in);
}

ImportanceSampler* make_fixed_epsilon_enki_abc_is(RandomNumberGenerator* rng_in,
                                                  size_t* seed_in,
                                                  Data* summary_data_in,
                                                  size_t number_of_particles_in,
                                                  IndependentProposalKernel* prior_in,
                                                  const std::vector<std::string> &summary_data_variables_in,
                                                  const std::string &scale_variable_in,
                                                  const std::string &epsilon_variable_in,
                                                  double min_epsilon,
                                                  size_t number_of_targets_in,
                                                  const arma::colvec &scale_in,
                                                  size_t enki_lag_in,
                                                  EnsembleShifter* shifter_in,
                                                  //double enki_annealing_desired_cess_in,
                                                  //size_t enki_number_of_bisections_in,
                                                  SimulateModelPtr simulate_model_in,
                                                  std::shared_ptr<Transform> summary_statistics_in,
                                                  bool enki_on_summary,
                                                  size_t number_of_abc_simulations_in,
                                                  double significance_level_in,
                                                  size_t estimator_type_in,
                                                  bool abc_parallel_in,
                                                  size_t abc_grain_size_in,
                                                  bool is_store_output_in,
                                                  bool is_parallel_in,
                                                  size_t is_grain_size_in,
                                                  const std::string &results_name_in)
{
  /*
   std::vector<double> enki_schedule;
   enki_schedule.push_back(arma::datum::inf);
   enki_schedule.push_back(epsilon_in);
   */
  
  std::vector<LikelihoodEstimator*> abc_likelihood;
  abc_likelihood.push_back(make_fixed_epsilon_enki_abc_likelihood(rng_in,
                                                                  seed_in,
                                                                  summary_data_in,
                                                                  enki_lag_in,
                                                                  shifter_in,
                                                                  //enki_annealing_desired_cess_in,
                                                                  //enki_number_of_bisections_in,
                                                                  scale_variable_in,
                                                                  epsilon_variable_in,
                                                                  min_epsilon,
                                                                  number_of_targets_in,
                                                                  scale_in,
                                                                  //enki_schedule,
                                                                  summary_data_variables_in,
                                                                  simulate_model_in,
                                                                  summary_statistics_in,
                                                                  enki_on_summary,
                                                                  number_of_abc_simulations_in,
                                                                  significance_level_in,
                                                                  estimator_type_in,
                                                                  abc_parallel_in,
                                                                  abc_grain_size_in));
  
  return new ImportanceSampler(rng_in,
                               seed_in,
                               summary_data_in,
                               Parameters(),
                               number_of_particles_in,
                               epsilon_variable_in,
                               min_epsilon,
                               abc_likelihood,
                               prior_in,
                               false,
                               is_store_output_in,
                               true,
                               true,
                               false,
                               is_parallel_in,
                               is_grain_size_in,
                               results_name_in);
}

SMCMCMCMove* make_fixed_epsilon_enki_abc_mcmc(RandomNumberGenerator* rng_in,
                                              size_t* seed_in,
                                              Data* data_in,
                                              MCMC* mcmc_in,
                                              const std::vector<Parameters> &initial_points_in,
                                              IndependentProposalKernel* prior_in,
                                              const std::vector<std::string> &data_variables_in,
                                              const std::string &scale_variable_in,
                                              const std::string &epsilon_variable_in,
                                              double min_epsilon,
                                              size_t number_of_targets_in,
                                              const arma::colvec &scale_in,
                                              size_t enki_lag_in,
                                              EnsembleShifter* shifter_in,
                                              //double enki_annealing_desired_cess_in,
                                              //size_t enki_number_of_bisections_in,
                                              SimulateModelPtr simulate_model_in,
                                              size_t number_of_abc_simulations_in,
                                              bool transform_proposed_particles,
                                              double significance_level_in,
                                              size_t estimator_type_in,
                                              bool abc_parallel_in,
                                              size_t abc_grain_size_in,
                                              bool mcmc_parallel_in,
                                              size_t mcmc_grain_size_in,
                                              const std::string &results_name_in)
{
  /*
   std::vector<double> enki_schedule;
   enki_schedule.push_back(arma::datum::inf);
   enki_schedule.push_back(epsilon_in);
   */
  
  std::vector<LikelihoodEstimator*> abc_likelihood_and_prior;
  abc_likelihood_and_prior.push_back(new ExactLikelihoodEstimator(rng_in,
                                                                  seed_in,
                                                                  data_in,
                                                                  prior_in,
                                                                  true));
  abc_likelihood_and_prior.push_back(make_fixed_epsilon_enki_abc_likelihood(rng_in,
                                                                            seed_in,
                                                                            data_in,
                                                                            enki_lag_in,
                                                                            shifter_in,
                                                                            //enki_annealing_desired_cess_in,
                                                                            //enki_number_of_bisections_in,
                                                                            scale_variable_in,
                                                                            epsilon_variable_in,
                                                                            min_epsilon,
                                                                            number_of_targets_in,
                                                                            scale_in,
                                                                            //enki_schedule,
                                                                            data_variables_in,
                                                                            simulate_model_in,
                                                                            number_of_abc_simulations_in,
                                                                            significance_level_in,
                                                                            estimator_type_in,
                                                                            abc_parallel_in,
                                                                            abc_grain_size_in));
  
  
  arma::colvec log_probabilities_of_initial_values(initial_points_in.size(),arma::fill::zeros);
  
  std::vector<size_t> indices;
  indices.reserve(abc_likelihood_and_prior.size());
  for (size_t i=0; i<abc_likelihood_and_prior.size(); ++i)
    indices.push_back(i);
  Index* an_index = new VectorIndex(indices);
  
  return new SMCMCMCMove(rng_in,
                         seed_in,
                         data_in,
                         Parameters(),
                         2,
                         2,
                         mcmc_in,
                         abc_likelihood_and_prior,
                         initial_points_in,
                         log_probabilities_of_initial_values,
                         an_index,
                         an_index->duplicate(),
                         transform_proposed_particles,
                         mcmc_parallel_in,
                         mcmc_grain_size_in,
                         results_name_in);
}

SMCMCMCMove* make_enkf_mcmc(RandomNumberGenerator* rng_in,
                            size_t* seed_in,
                            Data* data_in,
                            MCMC* mcmc_in,
                            const std::vector<Parameters> &initial_points_in,
                            IndependentProposalKernel* prior_in,
                            size_t enk_lag_in,
                            const std::string &index_name_in,
                            const std::string &time_name_in,
                            const std::string &time_diff_name_in,
                            size_t first_index_in,
                            size_t last_index_in,
                            size_t predictions_per_update_in,
                            double update_time_step_in,
                            double initial_time_in,
                            size_t number_of_ensemble_members_in,
                            EnsembleShifter* shifter_in,
                            IndependentProposalKernel* dynamic_prior_in,
                            ProposalKernel* transition_model_in,
                            const std::vector<MeasurementCovarianceEstimator*> &measurement_covariance_estimators_in,
                            bool transform_proposed_particles,
                            bool enkf_parallel_in,
                            size_t enkf_grain_size_in,
                            bool mcmc_parallel_in,
                            size_t mcmc_grain_size_in,
                            const std::string &results_name_in)
{
  std::vector<LikelihoodEstimator*> abc_likelihood_and_prior;
  abc_likelihood_and_prior.push_back(new ExactLikelihoodEstimator(rng_in,
                                                                  seed_in,
                                                                  data_in,
                                                                  prior_in,
                                                                  true));
  abc_likelihood_and_prior.push_back(make_enkf_likelihood(rng_in,
                                                          seed_in,
                                                          data_in,
                                                          enk_lag_in,
                                                          index_name_in,
                                                          time_name_in,
                                                          time_diff_name_in,
                                                          first_index_in,
                                                          last_index_in,
                                                          predictions_per_update_in,
                                                          update_time_step_in,
                                                          initial_time_in,
                                                          number_of_ensemble_members_in,
                                                          shifter_in,
                                                          dynamic_prior_in,
                                                          transition_model_in,
                                                          measurement_covariance_estimators_in,
                                                          enkf_parallel_in,
                                                          enkf_grain_size_in));
  
  
  arma::colvec log_probabilities_of_initial_values(initial_points_in.size(),arma::fill::zeros);
  
  std::vector<size_t> indices;
  indices.reserve(abc_likelihood_and_prior.size());
  for (size_t i=0; i<abc_likelihood_and_prior.size(); ++i)
    indices.push_back(i);
  Index* an_index = new VectorIndex(indices);
  
  return new SMCMCMCMove(rng_in,
                         seed_in,
                         data_in,
                         Parameters(),
                         2,
                         2,
                         mcmc_in,
                         abc_likelihood_and_prior,
                         initial_points_in,
                         log_probabilities_of_initial_values,
                         an_index,
                         an_index->duplicate(),
                         transform_proposed_particles,
                         mcmc_parallel_in,
                         mcmc_grain_size_in,
                         results_name_in);
}

ImportanceSampler* make_fixed_epsilon_gaussian_abc_is(RandomNumberGenerator* rng_in,
                                                      size_t* seed_in,
                                                      Data* summary_data_in,
                                                      size_t number_of_particles_in,
                                                      IndependentProposalKernel* prior_in,
                                                      const std::vector<std::string> &summary_data_variables_in,
                                                      const std::string &scale_variable_in,
                                                      const std::string &epsilon_variable_in,
                                                      double epsilon_in,
                                                      SimulateModelPtr simulate_model_in,
                                                      std::shared_ptr<Transform> summary_statistics_in,
                                                      size_t number_of_abc_simulations_in,
                                                      bool abc_store_output_in,
                                                      bool abc_parallel_in,
                                                      size_t abc_grain_size_in,
                                                      bool is_store_output_in,
                                                      bool is_parallel_in,
                                                      size_t is_grain_size_in,
                                                      const std::string &results_name_in)
{
  std::vector<LikelihoodEstimator*> abc_likelihood;
  abc_likelihood.push_back(make_fixed_epsilon_gaussian_abc_likelihood(rng_in,
                                                                      seed_in,
                                                                      summary_data_in,
                                                                      summary_data_variables_in,
                                                                      scale_variable_in,
                                                                      epsilon_variable_in,
                                                                      epsilon_in,
                                                                      simulate_model_in,
                                                                      summary_statistics_in,
                                                                      number_of_abc_simulations_in,
                                                                      abc_store_output_in,
                                                                      abc_parallel_in,
                                                                      abc_grain_size_in));
  
  return new ImportanceSampler(rng_in,
                               seed_in,
                               summary_data_in,
                               Parameters(),
                               number_of_particles_in,
                               epsilon_variable_in,
                               epsilon_in,
                               abc_likelihood,
                               prior_in,
                               false,
                               is_store_output_in,
                               true,
                               true,
                               false,
                               is_parallel_in,
                               is_grain_size_in,
                               results_name_in);
}

SMCMCMCMove* make_fixed_epsilon_lp_uniform_abc_mcmc(RandomNumberGenerator* rng_in,
                                                    size_t* seed_in,
                                                    Data* summary_data_in,
                                                    MCMC* mcmc_in,
                                                    const std::vector<Parameters> &initial_points_in,
                                                    DistributionFactor* prior_in,
                                                    double p_in,
                                                    const std::vector<std::string> &summary_data_variables_in,
                                                    const std::string &scale_variable_in,
                                                    const std::string &epsilon_variable_in,
                                                    double epsilon_in,
                                                    SimulateModelPtr simulate_model_in,
                                                    std::shared_ptr<Transform> summary_statistics_in,
                                                    size_t number_of_abc_simulations_in,
                                                    bool transform_proposed_particles,
                                                    bool abc_store_output_in,
                                                    bool abc_parallel_in,
                                                    size_t abc_grain_size_in,
                                                    bool mcmc_parallel_in,
                                                    size_t mcmc_grain_size_in,
                                                    const std::string &results_name_in)
{
  std::vector<LikelihoodEstimator*> abc_likelihood_and_prior;
  abc_likelihood_and_prior.push_back(new ExactLikelihoodEstimator(rng_in,
                                                                  seed_in,
                                                                  summary_data_in,
                                                                  prior_in,
                                                                  true));
  abc_likelihood_and_prior.push_back(make_fixed_epsilon_lp_uniform_abc_likelihood(rng_in,
                                                                                  seed_in,
                                                                                  summary_data_in,
                                                                                  p_in,
                                                                                  summary_data_variables_in,
                                                                                  scale_variable_in,
                                                                                  epsilon_variable_in,
                                                                                  epsilon_in,
                                                                                  simulate_model_in,
                                                                                  summary_statistics_in,
                                                                                  number_of_abc_simulations_in,
                                                                                  abc_store_output_in,
                                                                                  abc_parallel_in,
                                                                                  abc_grain_size_in));
  
  arma::colvec log_probabilities_of_initial_values(initial_points_in.size(),arma::fill::zeros);
  
  std::vector<size_t> indices;
  indices.reserve(abc_likelihood_and_prior.size());
  for (size_t i=0; i<abc_likelihood_and_prior.size(); ++i)
    indices.push_back(i);
  Index* an_index = new VectorIndex(indices);
  
  return new SMCMCMCMove(rng_in,
                         seed_in,
                         summary_data_in,
                         Parameters(),
                         2,
                         2,
                         mcmc_in,
                         abc_likelihood_and_prior,
                         initial_points_in,
                         log_probabilities_of_initial_values,
                         an_index,
                         an_index->duplicate(),
                         transform_proposed_particles,
                         mcmc_parallel_in,
                         mcmc_grain_size_in,
                         results_name_in);
}

SMCMCMCMove* make_fixed_epsilon_gaussian_abc_mcmc(RandomNumberGenerator* rng_in,
                                                  size_t* seed_in,
                                                  Data* data_in,
                                                  MCMC* mcmc_in,
                                                  const std::vector<Parameters> &initial_points_in,
                                                  IndependentProposalKernel* prior_in,
                                                  const std::vector<std::string> &data_variables_in,
                                                  const std::string &scale_variable_in,
                                                  const std::string &epsilon_variable_in,
                                                  double epsilon_in,
                                                  SimulateModelPtr simulate_model_in,
                                                  size_t number_of_abc_simulations_in,
                                                  bool transform_proposed_particles,
                                                  bool abc_store_output_in,
                                                  bool abc_parallel_in,
                                                  size_t abc_grain_size_in,
                                                  bool mcmc_parallel_in,
                                                  size_t mcmc_grain_size_in,
                                                  const std::string &results_name_in)
{
  std::vector<LikelihoodEstimator*> abc_likelihood_and_prior;
  abc_likelihood_and_prior.push_back(new ExactLikelihoodEstimator(rng_in,
                                                                  seed_in,
                                                                  data_in,
                                                                  prior_in,
                                                                  true));
  abc_likelihood_and_prior.push_back(make_fixed_epsilon_gaussian_abc_likelihood(rng_in,
                                                                                seed_in,
                                                                                data_in,
                                                                                data_variables_in,
                                                                                scale_variable_in,
                                                                                epsilon_variable_in,
                                                                                epsilon_in,
                                                                                simulate_model_in,
                                                                                number_of_abc_simulations_in,
                                                                                abc_store_output_in,
                                                                                abc_parallel_in,
                                                                                abc_grain_size_in));
  
  arma::colvec log_probabilities_of_initial_values(initial_points_in.size(),arma::fill::zeros);
  
  std::vector<size_t> indices;
  indices.reserve(abc_likelihood_and_prior.size());
  for (size_t i=0; i<abc_likelihood_and_prior.size(); ++i)
    indices.push_back(i);
  Index* an_index = new VectorIndex(indices);
  
  return new SMCMCMCMove(rng_in,
                         seed_in,
                         data_in,
                         Parameters(),
                         2,
                         2,
                         mcmc_in,
                         abc_likelihood_and_prior,
                         initial_points_in,
                         log_probabilities_of_initial_values,
                         an_index,
                         an_index->duplicate(),
                         transform_proposed_particles,
                         mcmc_parallel_in,
                         mcmc_grain_size_in,
                         results_name_in);
}

SMCMCMCMove* make_fixed_epsilon_gaussian_abc_mcmc(RandomNumberGenerator* rng_in,
                                                  size_t* seed_in,
                                                  Data* summary_data_in,
                                                  MCMC* mcmc_in,
                                                  const std::vector<Parameters> &initial_points_in,
                                                  IndependentProposalKernel* prior_in,
                                                  const std::vector<std::string> &summary_data_variables_in,
                                                  const std::string &scale_variable_in,
                                                  const std::string &epsilon_variable_in,
                                                  double epsilon_in,
                                                  SimulateModelPtr simulate_model_in,
                                                  std::shared_ptr<Transform> summary_statistics_in,
                                                  size_t number_of_abc_simulations_in,
                                                  bool transform_proposed_particles,
                                                  bool abc_store_output_in,
                                                  bool abc_parallel_in,
                                                  size_t abc_grain_size_in,
                                                  bool mcmc_parallel_in,
                                                  size_t mcmc_grain_size_in,
                                                  const std::string &results_name_in)
{
  std::vector<LikelihoodEstimator*> abc_likelihood_and_prior;
  abc_likelihood_and_prior.push_back(new ExactLikelihoodEstimator(rng_in,
                                                                  seed_in,
                                                                  summary_data_in,
                                                                  prior_in,
                                                                  true));
  abc_likelihood_and_prior.push_back(make_fixed_epsilon_gaussian_abc_likelihood(rng_in,
                                                                                seed_in,
                                                                                summary_data_in,
                                                                                summary_data_variables_in,
                                                                                scale_variable_in,
                                                                                epsilon_variable_in,
                                                                                epsilon_in,
                                                                                simulate_model_in,
                                                                                summary_statistics_in,
                                                                                number_of_abc_simulations_in,
                                                                                abc_store_output_in,
                                                                                abc_parallel_in,
                                                                                abc_grain_size_in));
  
  arma::colvec log_probabilities_of_initial_values(initial_points_in.size(),arma::fill::zeros);
  
  std::vector<size_t> indices;
  indices.reserve(abc_likelihood_and_prior.size());
  for (size_t i=0; i<abc_likelihood_and_prior.size(); ++i)
    indices.push_back(i);
  Index* an_index = new VectorIndex(indices);
  
  return new SMCMCMCMove(rng_in,
                         seed_in,
                         summary_data_in,
                         Parameters(),
                         2,
                         2,
                         mcmc_in,
                         abc_likelihood_and_prior,
                         initial_points_in,
                         log_probabilities_of_initial_values,
                         an_index,
                         an_index->duplicate(),
                         transform_proposed_particles,
                         mcmc_parallel_in,
                         mcmc_grain_size_in,
                         results_name_in);
}

SMCMCMCMove* make_ess_epsilon_lp_uniform_abc_is(RandomNumberGenerator* rng_in,
                                                size_t* seed_in,
                                                Data* data_in,
                                                size_t number_of_particles_in,
                                                IndependentProposalKernel* prior_in,
                                                double p_in,
                                                const std::vector<std::string> &data_variables_in,
                                                const std::string &scale_variable_in,
                                                const std::string &epsilon_variable_in,
                                                double target_ess_in,
                                                size_t number_of_bisections_in,
                                                SimulateModelPtr simulate_model_in,
                                                size_t number_of_abc_simulations_in,
                                                bool transform_proposed_particles,
                                                bool abc_store_output_in,
                                                bool abc_parallel_in,
                                                size_t abc_grain_size_in,
                                                bool is_store_output_in,
                                                bool is_parallel_in,
                                                size_t is_grain_size_in,
                                                const std::string &results_name_in)
{
  std::vector<LikelihoodEstimator*> abc_likelihood;
  /*
   abc_likelihood_and_prior.push_back(new ExactLikelihoodEstimator(rng_in,
   seed_in,
   data_in,
   prior_in,
   true));
   */
  abc_likelihood.push_back(make_varying_epsilon_lp_uniform_abc_likelihood(rng_in,
                                                                          seed_in,
                                                                          data_in,
                                                                          p_in,
                                                                          data_variables_in,
                                                                          scale_variable_in,
                                                                          epsilon_variable_in,
                                                                          simulate_model_in,
                                                                          number_of_abc_simulations_in,
                                                                          abc_store_output_in,
                                                                          abc_parallel_in,
                                                                          abc_grain_size_in));
  
  std::vector<double> sequence_values;
  sequence_values.push_back(arma::datum::inf);
  sequence_values.push_back(0.0);
  
  return new SMCMCMCMove(rng_in,
                         seed_in,
                         data_in,
                         Parameters(),
                         number_of_particles_in,
                         int(is_store_output_in),
                         0,
                         new MetropolisMCMC(),
                         0,
                         target_ess_in,
                         number_of_bisections_in,
                         new AlwaysSMCTermination(),
                         epsilon_variable_in,
                         sequence_values,
                         abc_likelihood,
                         prior_in,
                         false,
                         true,
                         true,
                         false,
                         transform_proposed_particles,
                         is_parallel_in,
                         is_grain_size_in,
                         results_name_in);
}

SMCMCMCMove* make_ess_epsilon_lp_uniform_abc_is(RandomNumberGenerator* rng_in,
                                                size_t* seed_in,
                                                Data* summary_data_in,
                                                size_t number_of_particles_in,
                                                IndependentProposalKernel* prior_in,
                                                double p_in,
                                                const std::vector<std::string> &summary_data_variables_in,
                                                const std::string &scale_variable_in,
                                                const std::string &epsilon_variable_in,
                                                double target_ess_in,
                                                size_t number_of_bisections_in,
                                                SimulateModelPtr simulate_model_in,
                                                std::shared_ptr<Transform> summary_statistics_in,
                                                size_t number_of_abc_simulations_in,
                                                bool transform_proposed_particles,
                                                bool abc_store_output_in,
                                                bool abc_parallel_in,
                                                size_t abc_grain_size_in,
                                                bool is_store_output_in,
                                                bool is_parallel_in,
                                                size_t is_grain_size_in,
                                                const std::string &results_name_in)
{
  std::vector<LikelihoodEstimator*> abc_likelihood;
  /*
   abc_likelihood_and_prior.push_back(new ExactLikelihoodEstimator(rng_in,
   seed_in,
   summary_data_in,
   prior_in,
   true));
   */
  abc_likelihood.push_back(make_varying_epsilon_lp_uniform_abc_likelihood(rng_in,
                                                                          seed_in,
                                                                          summary_data_in,
                                                                          p_in,
                                                                          summary_data_variables_in,
                                                                          scale_variable_in,
                                                                          epsilon_variable_in,
                                                                          simulate_model_in,
                                                                          summary_statistics_in,
                                                                          number_of_abc_simulations_in,
                                                                          abc_store_output_in,
                                                                          abc_parallel_in,
                                                                          abc_grain_size_in));
  
  std::vector<double> sequence_values;
  sequence_values.push_back(arma::datum::inf);
  sequence_values.push_back(0.0);
  
  return new SMCMCMCMove(rng_in,
                         seed_in,
                         summary_data_in,
                         Parameters(),
                         number_of_particles_in,
                         int(is_store_output_in),
                         0,
                         new MetropolisMCMC(),
                         0,
                         target_ess_in,
                         number_of_bisections_in,
                         new AlwaysSMCTermination(),
                         epsilon_variable_in,
                         sequence_values,
                         abc_likelihood,
                         prior_in,
                         false,
                         true,
                         true,
                         false,
                         transform_proposed_particles,
                         is_parallel_in,
                         is_grain_size_in,
                         results_name_in);
}

SMCMCMCMove* make_ess_epsilon_gaussian_abc_is(RandomNumberGenerator* rng_in,
                                              size_t* seed_in,
                                              Data* data_in,
                                              size_t number_of_particles_in,
                                              IndependentProposalKernel* prior_in,
                                              const std::vector<std::string> &data_variables_in,
                                              const std::string &scale_variable_in,
                                              const std::string &epsilon_variable_in,
                                              double target_ess_in,
                                              size_t number_of_bisections_in,
                                              SimulateModelPtr simulate_model_in,
                                              size_t number_of_abc_simulations_in,
                                              bool transform_proposed_particles,
                                              bool abc_store_output_in,
                                              bool abc_parallel_in,
                                              size_t abc_grain_size_in,
                                              bool is_store_output_in,
                                              bool is_parallel_in,
                                              size_t is_grain_size_in,
                                              const std::string &results_name_in)
{
  std::vector<LikelihoodEstimator*> abc_likelihood;
  /*
   abc_likelihood_and_prior.push_back(new ExactLikelihoodEstimator(rng_in,
   seed_in,
   data_in,
   prior_in,
   true));
   */
  abc_likelihood.push_back(make_varying_epsilon_gaussian_abc_likelihood(rng_in,
                                                                        seed_in,
                                                                        data_in,
                                                                        data_variables_in,
                                                                        scale_variable_in,
                                                                        epsilon_variable_in,
                                                                        simulate_model_in,
                                                                        number_of_abc_simulations_in,
                                                                        abc_store_output_in,
                                                                        abc_parallel_in,
                                                                        abc_grain_size_in));
  
  std::vector<double> sequence_values;
  sequence_values.push_back(arma::datum::inf);
  sequence_values.push_back(0.0);
  
  return new SMCMCMCMove(rng_in,
                         seed_in,
                         data_in,
                         Parameters(),
                         number_of_particles_in,
                         int(is_store_output_in),
                         0,
                         new MetropolisMCMC(),
                         0,
                         target_ess_in,
                         number_of_bisections_in,
                         new AlwaysSMCTermination(),
                         epsilon_variable_in,
                         sequence_values,
                         abc_likelihood,
                         prior_in,
                         false,
                         true,
                         true,
                         false,
                         transform_proposed_particles,
                         is_parallel_in,
                         is_grain_size_in,
                         results_name_in);
}

SMCMCMCMove* make_ess_epsilon_gaussian_abc_is(RandomNumberGenerator* rng_in,
                                              size_t* seed_in,
                                              Data* summary_data_in,
                                              size_t number_of_particles_in,
                                              IndependentProposalKernel* prior_in,
                                              const std::vector<std::string> &summary_data_variables_in,
                                              const std::string &scale_variable_in,
                                              const std::string &epsilon_variable_in,
                                              double target_ess_in,
                                              size_t number_of_bisections_in,
                                              SimulateModelPtr simulate_model_in,
                                              std::shared_ptr<Transform> summary_statistics_in,
                                              size_t number_of_abc_simulations_in,
                                              bool transform_proposed_particles,
                                              bool abc_store_output_in,
                                              bool abc_parallel_in,
                                              size_t abc_grain_size_in,
                                              bool is_store_output_in,
                                              bool is_parallel_in,
                                              size_t is_grain_size_in,
                                              const std::string &results_name_in)
{
  std::vector<LikelihoodEstimator*> abc_likelihood;
  /*
   abc_likelihood_and_prior.push_back(new ExactLikelihoodEstimator(rng_in,
   seed_in,
   summary_data_in,
   prior_in,
   true));
   */
  abc_likelihood.push_back(make_varying_epsilon_gaussian_abc_likelihood(rng_in,
                                                                        seed_in,
                                                                        summary_data_in,
                                                                        summary_data_variables_in,
                                                                        scale_variable_in,
                                                                        epsilon_variable_in,
                                                                        simulate_model_in,
                                                                        summary_statistics_in,
                                                                        number_of_abc_simulations_in,
                                                                        abc_store_output_in,
                                                                        abc_parallel_in,
                                                                        abc_grain_size_in));
  
  std::vector<double> sequence_values;
  sequence_values.push_back(arma::datum::inf);
  sequence_values.push_back(0.0);
  
  return new SMCMCMCMove(rng_in,
                         seed_in,
                         summary_data_in,
                         Parameters(),
                         number_of_particles_in,
                         int(is_store_output_in),
                         0,
                         new MetropolisMCMC(),
                         0,
                         target_ess_in,
                         number_of_bisections_in,
                         new AlwaysSMCTermination(),
                         epsilon_variable_in,
                         sequence_values,
                         abc_likelihood,
                         prior_in,
                         false,
                         true,
                         true,
                         false,
                         transform_proposed_particles,
                         is_parallel_in,
                         is_grain_size_in,
                         results_name_in);
}

SMCMCMCMove* make_cess_epsilon_lp_uniform_abc_smc(RandomNumberGenerator* rng_in,
                                                  size_t* seed_in,
                                                  Data* summary_data_in,
                                                  size_t number_of_particles_in,
                                                  MCMC* mcmc_in,
                                                  double resampling_desired_ess_in,
                                                  IndependentProposalKernel* prior_in,
                                                  double p_in,
                                                  const std::vector<std::string> &summary_data_variables_in,
                                                  const std::string &scale_variable_in,
                                                  const std::string &epsilon_variable_in,
                                                  const std::vector<double> &schedule_in,
                                                  double target_cess_in,
                                                  size_t number_of_bisections_in,
                                                  size_t number_of_iterations_for_similar_epsilon_for_termination,
                                                  double threshold_to_determine_similar_epsilon_for_termination,
                                                  SimulateModelPtr simulate_model_in,
                                                  std::shared_ptr<Transform> summary_statistics_in,
                                                  size_t number_of_abc_simulations_in,
                                                  bool transform_proposed_particles,
                                                  bool abc_store_output_in,
                                                  bool abc_parallel_in,
                                                  size_t abc_grain_size_in,
                                                  bool smc_parallel_in,
                                                  size_t smc_grain_size_in,
                                                  const std::string &results_name_in)
{
  std::vector<LikelihoodEstimator*> abc_likelihood_and_prior;
  abc_likelihood_and_prior.push_back(new ExactLikelihoodEstimator(rng_in,
                                                                  seed_in,
                                                                  summary_data_in,
                                                                  prior_in,
                                                                  true));
  abc_likelihood_and_prior.push_back(make_varying_epsilon_lp_uniform_abc_likelihood(rng_in,
                                                                                    seed_in,
                                                                                    summary_data_in,
                                                                                    p_in,
                                                                                    summary_data_variables_in,
                                                                                    scale_variable_in,
                                                                                    epsilon_variable_in,
                                                                                    simulate_model_in,
                                                                                    summary_statistics_in,
                                                                                    number_of_abc_simulations_in,
                                                                                    abc_store_output_in,
                                                                                    abc_parallel_in,
                                                                                    abc_grain_size_in));
  
  std::vector<double> sequence_values = schedule_in;
  
  std::cout << "need to read in lag if we ever use this function for something other than testing." << std::endl;
  
  return new SMCMCMCMove(rng_in,
                         seed_in,
                         summary_data_in,
                         Parameters(),
                         number_of_particles_in,
                         10000000000,
                         10000000000,
                         mcmc_in,
                         resampling_desired_ess_in,
                         target_cess_in,
                         number_of_bisections_in,
                         new StableSMCTermination(number_of_iterations_for_similar_epsilon_for_termination,
                                                  threshold_to_determine_similar_epsilon_for_termination),
                         epsilon_variable_in,
                         sequence_values,
                         abc_likelihood_and_prior,
                         prior_in,
                         true,
                         true,
                         true,
                         false,
                         transform_proposed_particles,
                         smc_parallel_in,
                         smc_grain_size_in,
                         results_name_in);
}

ImportanceSampler* make_sl_is(RandomNumberGenerator* rng_in,
                              size_t* seed_in,
                              Data* summary_data_in,
                              size_t number_of_particles_in,
                              IndependentProposalKernel* prior_in,
                              const std::vector<std::string> &summary_data_variables_in,
                              SimulateModelPtr simulate_model_in,
                              std::shared_ptr<Transform> summary_statistics_in,
                              size_t number_of_sl_simulations_in,
                              bool unbiased_in,
                              bool store_output_in,
                              bool sl_store_output_in,
                              bool sl_parallel_in,
                              size_t sl_grain_size_in,
                              bool is_store_output_in,
                              bool is_parallel_in,
                              size_t is_grain_size_in,
                              const std::string &results_name_in)
{
  std::vector<LikelihoodEstimator*> sl_likelihood;
  sl_likelihood.push_back(make_sl_likelihood(rng_in,
                                             seed_in,
                                             summary_data_in,
                                             summary_data_variables_in,
                                             simulate_model_in,
                                             summary_statistics_in,
                                             number_of_sl_simulations_in,
                                             unbiased_in,
                                             sl_store_output_in,
                                             sl_parallel_in,
                                             sl_grain_size_in));
  
  return new ImportanceSampler(rng_in,
                               seed_in,
                               summary_data_in,
                               Parameters(),
                               number_of_particles_in,
                               "",
                               sl_likelihood,
                               prior_in,
                               false,
                               is_store_output_in,
                               true,
                               true,
                               false,
                               is_parallel_in,
                               is_grain_size_in,
                               results_name_in);
}

SMCMCMCMove* make_sl_mcmc(RandomNumberGenerator* rng_in,
                          size_t* seed_in,
                          Data* summary_data_in,
                          MCMC* mcmc_in,
                          const std::vector<Parameters> &initial_points_in,
                          DistributionFactor* prior_in,
                          const std::vector<std::string> &summary_data_variables_in,
                          SimulateModelPtr simulate_model_in,
                          std::shared_ptr<Transform> summary_statistics_in,
                          size_t number_of_sl_simulations_in,
                          bool unbiased_in,
                          bool transform_proposed_particles,
                          bool sl_store_output_in,
                          bool sl_parallel_in,
                          size_t sl_grain_size_in,
                          bool mcmc_parallel_in,
                          size_t mcmc_grain_size_in,
                          const std::string &results_name_in)
{
  std::vector<LikelihoodEstimator*> sl_likelihood_and_prior;
  sl_likelihood_and_prior.push_back(new ExactLikelihoodEstimator(rng_in,
                                                                 seed_in,
                                                                 summary_data_in,
                                                                 prior_in,
                                                                 true));
  sl_likelihood_and_prior.push_back(make_sl_likelihood(rng_in,
                                                       seed_in,
                                                       summary_data_in,
                                                       summary_data_variables_in,
                                                       simulate_model_in,
                                                       summary_statistics_in,
                                                       number_of_sl_simulations_in,
                                                       unbiased_in,
                                                       sl_store_output_in,
                                                       sl_parallel_in,
                                                       sl_grain_size_in));
  
  arma::colvec log_probabilities_of_initial_values(initial_points_in.size(),arma::fill::zeros);
  
  std::vector<size_t> indices;
  indices.reserve(sl_likelihood_and_prior.size());
  for (size_t i=0; i<sl_likelihood_and_prior.size(); ++i)
    indices.push_back(i);
  Index* an_index = new VectorIndex(indices);
  
  return new SMCMCMCMove(rng_in,
                         seed_in,
                         summary_data_in,
                         Parameters(),
                         2,
                         2,
                         mcmc_in,
                         sl_likelihood_and_prior,
                         initial_points_in,
                         log_probabilities_of_initial_values,
                         an_index,
                         an_index->duplicate(),
                         transform_proposed_particles,
                         mcmc_parallel_in,
                         mcmc_grain_size_in,
                         results_name_in);
}

SMCMCMCMove* make_cess_annealing_sl_smc(RandomNumberGenerator* rng_in,
                                        size_t* seed_in,
                                        Data* summary_data_in,
                                        size_t number_of_particles_in,
                                        MCMC* mcmc_in,
                                        double resampling_desired_ess_in,
                                        IndependentProposalKernel* prior_in,
                                        const std::vector<std::string> &summary_data_variables_in,
                                        const std::string &annealing_variable_in,
                                        double annealing_desired_cess_in,
                                        size_t number_of_bisections_in,
                                        SimulateModelPtr simulate_model_in,
                                        std::shared_ptr<Transform> summary_statistics_in,
                                        size_t number_of_sl_simulations_in,
                                        bool unbiased_in,
                                        bool transform_proposed_particles,
                                        bool sl_store_output_in,
                                        bool sl_parallel_in,
                                        size_t sl_grain_size_in,
                                        bool smc_parallel_in,
                                        size_t smc_grain_size_in,
                                        const std::string &results_name_in)
{
  std::vector<LikelihoodEstimator*> annealed_sl_likelihood_and_prior;
  annealed_sl_likelihood_and_prior.push_back(new ExactLikelihoodEstimator(rng_in,
                                                                          seed_in,
                                                                          summary_data_in,
                                                                          prior_in,
                                                                          true));
  DensityLikelihoodEstimator* sl = make_sl_likelihood(rng_in,
                                                      seed_in,
                                                      summary_data_in,
                                                      summary_data_variables_in,
                                                      simulate_model_in,
                                                      summary_statistics_in,
                                                      number_of_sl_simulations_in,
                                                      unbiased_in,
                                                      sl_store_output_in,
                                                      sl_parallel_in,
                                                      sl_grain_size_in);
  annealed_sl_likelihood_and_prior.push_back(new AnnealedLikelihoodEstimator(rng_in,
                                                                             seed_in,
                                                                             summary_data_in,
                                                                             sl,
                                                                             annealing_power,
                                                                             annealing_variable_in,
                                                                             false));
  
  std::vector<double> temperatures;
  temperatures.push_back(0.0);
  temperatures.push_back(1.0);
  
  return new SMCMCMCMove(rng_in,
                         seed_in,
                         summary_data_in,
                         Parameters(),
                         number_of_particles_in,
                         2,
                         2,
                         mcmc_in,
                         resampling_desired_ess_in,
                         annealing_desired_cess_in,
                         number_of_bisections_in,
                         NULL,
                         annealing_variable_in,
                         temperatures,
                         annealed_sl_likelihood_and_prior,
                         prior_in,
                         true,
                         true,
                         true,
                         false,
                         transform_proposed_particles,
                         smc_parallel_in,
                         smc_grain_size_in,
                         results_name_in);
}

/*
 LikelihoodEstimator* make_likelihood_estimator(const List &model,
 const List &algorithm)
 {
 std::string likelihood_method = Rcpp::as<std::string>(algorithm["likelihood_method"]);
 // if (likelihood_method=="analytic")
 // {
 //   SEXP evaluate_log_likelihood_SEXP = model["evaluate_log_likelihood"];
 //   EvaluateLogLikelihoodPtr evaluate_log_likelihood = load_evaluate_log_likelihood(evaluate_log_likelihood_SEXP);
 //
 //   List observed_data = model["observed_data"];
 //
 //   return new ExactLikelihoodEstimator(observed_data, evaluate_log_likelihood);
 // }
 // else if (likelihood_method=="abc")
 // {
 //   SEXP simulate_model_SEXP = model["simulate_model"];
 //   SimulateModelPtr simulate_model = load_simulate_model(simulate_model_SEXP);
 //
 //   //SEXP get_data_from_simulation_SEXP = algorithm["get_data_from_simulation"];
 //   //GetDataFromSimulationPtr get_data_from_simulation = load_get_data_from_simulation(get_data_from_simulation_SEXP);
 //
 //   unsigned int number_of_likelihood_particles = algorithm["number_of_likelihood_particles"];
 //
 //   SEXP evaluate_log_abc_kernel_SEXP = algorithm["evaluate_log_abc_kernel"];
 //   EvaluateLogABCKernelPtr evaluate_log_abc_kernel = load_evaluate_log_abc_kernel(evaluate_log_abc_kernel_SEXP);
 //
 //   SEXP summary_statistics_SEXP = algorithm["summary_statistics"];
 //   std::shared_ptr<Transform> summary_statistics = load_summary_statistics(summary_statistics_SEXP);
 //
 //   double abc_tolerance = algorithm["abc_tolerance"];
 //
 //   double abc_desired_cess = algorithm["abc_desired_cess"];
 //
 //   arma::colvec summary_statistics_scaling = algorithm["summary_statistics_scaling"];
 //
 //   bool adapt_abc_tolerance_to_cess = algorithm["adapt_abc_tolerance_to_cess"];
 //
 //   bool adapt_summary_statistics_scaling = algorithm["adapt_summary_statistics_scaling"];
 //
 //   List observed_data = model["observed_data"];
 //
 //   return new ABCLikelihoodEstimator(observed_data,
 //                                     simulate_model,
 //                                     number_of_likelihood_particles,
 //                                     evaluate_log_abc_kernel,
 //                                     summary_statistics,
 //                                     abc_tolerance,
 //                                     abc_desired_cess,
 //                                     summary_statistics_scaling,
 //                                     adapt_abc_tolerance_to_cess,
 //                                     adapt_summary_statistics_scaling);
 // }
 
 Rcpp::stop("likelihood_method not set to a valid option.");
 }
 */
}
