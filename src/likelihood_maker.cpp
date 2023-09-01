#include "ilike_header.h"
#include "exact_likelihood_estimator.h"
#include "abc_likelihood_estimator.h"
#include "lp_uniform_abc_kernel_factor.h"
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

ABCLikelihoodEstimator* make_lp_uniform_abc_kernel(RandomNumberGenerator* rng_in,
                                                   size_t* seed_in,
                                                   Data* data_in,
                                                   double p_in,
                                                   const std::vector<std::string> &data_variables_in,
                                                   const std::string &epsilon_variable_in,
                                                   bool fixed_epsilon_in)
{
  return new ABCLikelihoodEstimator(rng_in,
                                    seed_in,
                                    data_in,
                                    new LpUniformABCKernelFactor(p_in,
                                                                 data_variables_in,
                                                                 epsilon_variable_in,
                                                                 data_in),
                                    fixed_epsilon_in);
}

ImportanceSampler* make_fixed_epsilon_lp_uniform_abc_likelihood(RandomNumberGenerator* rng_in,
                                                                size_t* seed_in,
                                                                Data* data_in,
                                                                double p_in,
                                                                const std::vector<std::string> &data_variables_in,
                                                                const std::string &epsilon_variable_in,
                                                                double epsilon_in,
                                                                SimulateModelPtr simulate_model_in,
                                                                size_t number_of_abc_simulations_in,
                                                                bool parallel_in,
                                                                size_t grain_size_in)
{
  std::vector<LikelihoodEstimator*> abc_kernel;
  abc_kernel.push_back(make_lp_uniform_abc_kernel(rng_in,
                                                  seed_in,
                                                  data_in,
                                                  p_in,
                                                  data_variables_in,
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
                                                                const std::string &epsilon_variable_in,
                                                                double epsilon_in,
                                                                SimulateModelPtr simulate_model_in,
                                                                SummaryStatisticsPtr summary_statistics_in,
                                                                size_t number_of_abc_simulations_in,
                                                                bool parallel_in,
                                                                size_t grain_size_in)
{
  std::vector<LikelihoodEstimator*> abc_kernel;
  abc_kernel.push_back(make_lp_uniform_abc_kernel(rng_in,
                                                  seed_in,
                                                  summary_data_in,
                                                  p_in,
                                                  summary_data_variables_in,
                                                  epsilon_variable_in,
                                                  true));
  
  IndependentProposalKernel* model_simulator = new CustomIndependentProposalKernel(simulate_model_in);
  IndependentProposalKernel* summary_model_simulator = new TransformedIndependentProposalKernel(model_simulator,
                                                                                                Transform(summary_statistics_in),
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
                                                                  const std::string &epsilon_variable_in,
                                                                  SimulateModelPtr simulate_model_in,
                                                                  size_t number_of_abc_simulations_in,
                                                                  bool parallel_in,
                                                                  size_t grain_size_in)
{
  std::vector<LikelihoodEstimator*> abc_kernel;
  abc_kernel.push_back(make_lp_uniform_abc_kernel(rng_in,
                                                  seed_in,
                                                  data_in,
                                                  p_in,
                                                  data_variables_in,
                                                  epsilon_variable_in,
                                                  false));
  
  IndependentProposalKernel* model_simulator = new CustomIndependentProposalKernel(simulate_model_in);
  
  return new ImportanceSampler(rng_in,
                               seed_in,
                               data_in,
                               Parameters(),
                               number_of_abc_simulations_in,
                               abc_kernel,
                               model_simulator,
                               false,
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
                                                                  const std::string &epsilon_variable_in,
                                                                  SimulateModelPtr simulate_model_in,
                                                                  SummaryStatisticsPtr summary_statistics_in,
                                                                  size_t number_of_abc_simulations_in,
                                                                  bool parallel_in,
                                                                  size_t grain_size_in)
{
  std::vector<LikelihoodEstimator*> abc_kernel;
  abc_kernel.push_back(make_lp_uniform_abc_kernel(rng_in,
                                                  seed_in,
                                                  summary_data_in,
                                                  p_in,
                                                  summary_data_variables_in,
                                                  epsilon_variable_in,
                                                  false));
  
  IndependentProposalKernel* model_simulator = new CustomIndependentProposalKernel(simulate_model_in);
  IndependentProposalKernel* summary_model_simulator = new TransformedIndependentProposalKernel(model_simulator,
                                                                                                Transform(summary_statistics_in));
  
  return new ImportanceSampler(rng_in,
                               seed_in,
                               summary_data_in,
                               Parameters(),
                               number_of_abc_simulations_in,
                               abc_kernel,
                               summary_model_simulator,
                               false,
                               false,
                               true,
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
                                        parallel_in,
                                        grain_size_in);
}

DensityLikelihoodEstimator* make_sl_likelihood(RandomNumberGenerator* rng_in,
                                               size_t* seed_in,
                                               Data* summary_data_in,
                                               const std::vector<std::string> &summary_data_variables_in,
                                               SimulateModelPtr simulate_model_in,
                                               SummaryStatisticsPtr summary_statistics_in,
                                               size_t number_of_sl_simulations_in,
                                               bool unbiased_in,
                                               bool parallel_in,
                                               size_t grain_size_in)
{
  IndependentProposalKernel* model_simulator = new CustomIndependentProposalKernel(simulate_model_in);
  IndependentProposalKernel* summary_model_simulator = new TransformedIndependentProposalKernel(model_simulator,
                                                                                                Transform(summary_statistics_in));
  
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
                                                        const std::string &epsilon_variable_in,
                                                        double epsilon_in,
                                                        SimulateModelPtr simulate_model_in,
                                                        size_t number_of_abc_simulations_in,
                                                        bool abc_parallel_in,
                                                        size_t abc_grain_size_in,
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
                                                                        epsilon_variable_in,
                                                                        epsilon_in,
                                                                        simulate_model_in,
                                                                        number_of_abc_simulations_in,
                                                                        abc_parallel_in,
                                                                        abc_grain_size_in));
  
  return new ImportanceSampler(rng_in,
                               seed_in,
                               data_in,
                               Parameters(),
                               number_of_particles_in,
                               abc_likelihood,
                               prior_in,
                               false,
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
                                                        const std::string &epsilon_variable_in,
                                                        double epsilon_in,
                                                        SimulateModelPtr simulate_model_in,
                                                        SummaryStatisticsPtr summary_statistics_in,
                                                        size_t number_of_abc_simulations_in,
                                                        bool abc_parallel_in,
                                                        size_t abc_grain_size_in,
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
                                                                        epsilon_variable_in,
                                                                        epsilon_in,
                                                                        simulate_model_in,
                                                                        summary_statistics_in,
                                                                        number_of_abc_simulations_in,
                                                                        abc_parallel_in,
                                                                        abc_grain_size_in));
  
  return new ImportanceSampler(rng_in,
                               seed_in,
                               summary_data_in,
                               Parameters(),
                               number_of_particles_in,
                               abc_likelihood,
                               prior_in,
                               false,
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
                                                    const std::string &epsilon_variable_in,
                                                    double epsilon_in,
                                                    SimulateModelPtr simulate_model_in,
                                                    SummaryStatisticsPtr summary_statistics_in,
                                                    size_t number_of_abc_simulations_in,
                                                    bool transform_proposed_particles,
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
                                                                                  epsilon_variable_in,
                                                                                  epsilon_in,
                                                                                  simulate_model_in,
                                                                                  summary_statistics_in,
                                                                                  number_of_abc_simulations_in,
                                                                                  abc_parallel_in,
                                                                                  abc_grain_size_in));
  
  arma::colvec log_probabilities_of_initial_values(initial_points_in.size(),arma::fill::zeros);

  return new SMCMCMCMove(rng_in,
                         seed_in,
                         summary_data_in,
                         Parameters(),
                         0,
                         0,
                         mcmc_in,
                         abc_likelihood_and_prior,
                         initial_points_in,
                         log_probabilities_of_initial_values,
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
                                                 const std::string &epsilon_variable_in,
                                                 double target_ess_in,
                                                 size_t number_of_bisections_in,
                                                 SimulateModelPtr simulate_model_in,
                                                 size_t number_of_abc_simulations_in,
                                                 bool transform_proposed_particles,
                                                 bool abc_parallel_in,
                                                 size_t abc_grain_size_in,
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
                                                                          epsilon_variable_in,
                                                                          simulate_model_in,
                                                                          number_of_abc_simulations_in,
                                                                          abc_parallel_in,
                                                                          abc_grain_size_in));
  
  std::vector<double> sequence_values;
  sequence_values.push_back(0.0);
  
  return new SMCMCMCMove(rng_in,
                         seed_in,
                         data_in,
                         Parameters(),
                         number_of_particles_in,
                         0,
                         0,
                         NULL,
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
                                                const std::string &epsilon_variable_in,
                                                double target_ess_in,
                                                size_t number_of_bisections_in,
                                                SimulateModelPtr simulate_model_in,
                                                SummaryStatisticsPtr summary_statistics_in,
                                                size_t number_of_abc_simulations_in,
                                                bool transform_proposed_particles,
                                                bool abc_parallel_in,
                                                size_t abc_grain_size_in,
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
                                                                          epsilon_variable_in,
                                                                          simulate_model_in,
                                                                          summary_statistics_in,
                                                                          number_of_abc_simulations_in,
                                                                          abc_parallel_in,
                                                                          abc_grain_size_in));
  
  std::vector<double> sequence_values;
  sequence_values.push_back(0.0);
  
  return new SMCMCMCMove(rng_in,
                         seed_in,
                         summary_data_in,
                         Parameters(),
                         number_of_particles_in,
                         0,
                         0,
                         NULL,
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
                                                  const std::string &epsilon_variable_in,
                                                  double target_cess_in,
                                                  size_t number_of_bisections_in,
                                                  size_t number_of_iterations_for_similar_epsilon_for_termination,
                                                  double threshold_to_determine_similar_epsilon_for_termination,
                                                  SimulateModelPtr simulate_model_in,
                                                  SummaryStatisticsPtr summary_statistics_in,
                                                  size_t number_of_abc_simulations_in,
                                                  bool transform_proposed_particles,
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
                                                                                    epsilon_variable_in,
                                                                                    simulate_model_in,
                                                                                    summary_statistics_in,
                                                                                    number_of_abc_simulations_in,
                                                                                    abc_parallel_in,
                                                                                    abc_grain_size_in));
  
  std::vector<double> sequence_values;
  sequence_values.push_back(0.0);
  
  return new SMCMCMCMove(rng_in,
                         seed_in,
                         summary_data_in,
                         Parameters(),
                         number_of_particles_in,
                         0,
                         0,
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
                              SummaryStatisticsPtr summary_statistics_in,
                              size_t number_of_sl_simulations_in,
                              bool unbiased_in,
                              bool sl_parallel_in,
                              size_t sl_grain_size_in,
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
                                             sl_parallel_in,
                                             sl_grain_size_in));
  
  return new ImportanceSampler(rng_in,
                               seed_in,
                               summary_data_in,
                               Parameters(),
                               number_of_particles_in,
                               sl_likelihood,
                               prior_in,
                               false,
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
                          SummaryStatisticsPtr summary_statistics_in,
                          size_t number_of_sl_simulations_in,
                          bool unbiased_in,
                          bool transform_proposed_particles,
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
                                                       sl_parallel_in,
                                                       sl_grain_size_in));
  
  arma::colvec log_probabilities_of_initial_values(initial_points_in.size(),arma::fill::zeros);
  
  return new SMCMCMCMove(rng_in,
                         seed_in,
                         summary_data_in,
                         Parameters(),
                         0,
                         0,
                         mcmc_in,
                         sl_likelihood_and_prior,
                         initial_points_in,
                         log_probabilities_of_initial_values,
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
                                        SummaryStatisticsPtr summary_statistics_in,
                                        size_t number_of_sl_simulations_in,
                                        bool unbiased_in,
                                        bool transform_proposed_particles,
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
                         0,
                         0,
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
  //   SummaryStatisticsPtr summary_statistics = load_summary_statistics(summary_statistics_SEXP);
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
