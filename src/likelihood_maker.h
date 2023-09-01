#ifndef LIKELIHOOD_MAKER
#define LIKELIHOOD_MAKER

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "likelihood_estimator.h"

class ABCLikelihoodEstimator;
class DensityLikelihoodEstimator;
class ImportanceSampler;
class SMCMCMCMove;


// ABC kernels

ABCLikelihoodEstimator* make_lp_uniform_abc_kernel(RandomNumberGenerator* rng_in,
                                                   size_t* seed_in,
                                                   Data* data_in,
                                                   double p_in,
                                                   const std::vector<std::string> &data_variables_in,
                                                   const std::string &epsilon_variable_in,
                                                   bool fixed_epsilon_in);


// ABC likelihoods

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
                                                                size_t grain_size_in);

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
                                                                size_t grain_size_in);

ImportanceSampler* make_varying_epsilon_lp_uniform_abc_likelihood(RandomNumberGenerator* rng_in,
                                                                  size_t* seed_in,
                                                                  Data* data_in,
                                                                  double p_in,
                                                                  const std::vector<std::string> &data_variables_in,
                                                                  const std::string &epsilon_variable_in,
                                                                  SimulateModelPtr simulate_model_in,
                                                                  size_t number_of_abc_simulations_in,
                                                                  bool parallel_in,
                                                                  size_t grain_size_in);

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
                                                                  size_t grain_size_in);


// Synthetic likelihoods.

DensityLikelihoodEstimator* make_sl_likelihood(RandomNumberGenerator* rng_in,
                                               size_t* seed_in,
                                               Data* data_in,
                                               const std::vector<std::string> &data_variables_in,
                                               SimulateModelPtr simulate_model_in,
                                               size_t number_of_sl_simulations_in,
                                               bool unbiased_in,
                                               bool parallel_in,
                                               size_t grain_size_in);

DensityLikelihoodEstimator* make_sl_likelihood(RandomNumberGenerator* rng_in,
                                               size_t* seed_in,
                                               Data* data_in,
                                               const std::vector<std::string> &data_variables_in,
                                               SimulateModelPtr simulate_model_in,
                                               SummaryStatisticsPtr summary_statistics_in,
                                               size_t number_of_sl_simulations_in,
                                               bool unbiased_in,
                                               bool parallel_in,
                                               size_t grain_size_in);

// ABC samplers

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
                                                        const std::string &results_name_in);

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
                                                        const std::string &results_name_in);

// abc-mcmc
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
                                                    bool abc_parallel_in,
                                                    size_t abc_grain_size_in,
                                                    bool mcmc_parallel_in,
                                                    size_t mcmc_grain_size_in,
                                                    const std::string &results_name_in);

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
                                                 bool abc_parallel_in,
                                                 size_t abc_grain_size_in,
                                                 bool is_parallel_in,
                                                 size_t is_grain_size_in,
                                                 const std::string &results_name_in);

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
                                                bool abc_parallel_in,
                                                size_t abc_grain_size_in,
                                                bool is_parallel_in,
                                                size_t is_grain_size_in,
                                                const std::string &results_name_in);

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
                                                  bool abc_parallel_in,
                                                  size_t abc_grain_size_in,
                                                  bool smc_parallel_in,
                                                  size_t smc_grain_size_in,
                                                  const std::string &results_name_in);

// SL samplers.

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
                              const std::string &results_name_in);

// sl-mcmc
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
                          bool sl_parallel_in,
                          size_t sl_grain_size_in,
                          bool mcmc_parallel_in,
                          size_t mcmc_grain_size_in,
                          const std::string &results_name_in);

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
                                         const std::string &results_name_in);


#endif
