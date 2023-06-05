#ifndef SMC_H
#define SMC_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include "model_and_algorithm.h"
#include "likelihood_estimator.h"
#include "particles.h"
#include "distributions.h"
#include "parameters.h"
#include "sequencer.h"
#include "ilike_header.h"

class SMCOutput;
class SMCWorker;
class SequentialSMCWorker;
class SMCCriterion;
class MCMC;
class MoveOutput;
class ParticleSimulator;

class SMC : public LikelihoodEstimator
{
public:

  SMC();
  SMC(RandomNumberGenerator* rng_in,
      size_t* seed_in,
      Data* data_in,
      const Parameters &algorithm_parameters_in,
      size_t number_of_particles_in,
      size_t lag_in,
      size_t lag_proposed_in,
      //SummaryStatisticsPtr summary_statistics_in,
      double resampling_desired_ess_in,
      //EvaluateLogDistributionPtr evaluate_log_proposal_in,
      bool proposal_is_evaluated_in,
      bool smcfixed_flag_in,
      bool sequencer_limit_is_fixed_in,
      const std::string &results_name_in);
  
  //SMC(RandomNumberGenerator* rng_in,
  //    size_t* seed_in,
  //    Data* data_in,
  //    size_t lag_in,
  //    size_t lag_proposed_in,
  //    double resampling_desired_ess_in,
  //    const std::vector<Parameters> &initial_values_in,
  //    const arma::colvec &log_probabilities_of_initial_values_in,
  //    bool smcfixed_flag_in,
  //    bool sequencer_limit_is_fixed_in);
  
  virtual ~SMC();

  SMC(const SMC &another);
  void operator=(const SMC &another);
  virtual SMC* smc_duplicate() const=0;

  SMCOutput* run();
  SMCOutput* run(const Parameters &conditioned_on_parameters);
  
  SMCOutput* run(const std::string &directory_name);
  SMCOutput* run(const std::string &directory_name,
                 const Parameters &conditioned_on_parameters);

  LikelihoodEstimatorOutput* initialise();
  LikelihoodEstimatorOutput* initialise(const Parameters &conditioned_on_parameters);
  
  void setup();
  void setup(const Parameters &parameters);
  void setup_variables();
  void setup_variables(const Parameters &parameters);
  
  void resample(SMCOutput* current_state);
  
  virtual MoveOutput* move(RandomNumberGenerator &rng,
                           Particle &particle)=0;
  //virtual void specific_move(RandomNumberGenerator &rng,
  //                           Particle* new_particle,
  //                           Particle &particle)=0;
  
  //virtual void weight_for_adapting_sequence(const Index* index,
  //                                          Particles &current_particles)=0;
  
  /*
  virtual MoveOutput* move(RandomNumberGenerator &rng,
                           Particle &particle,
                           const Parameters &conditioned_on_parameters)=0;
  */
  //virtual void specific_move(RandomNumberGenerator &rng,
  //                           Particle* new_particle,
  //                           Particle &particle,
  //                           const Parameters &conditioned_on_parameters)=0;
  
  virtual void weight_for_adapting_sequence(const Index* index,
                                            Particles &current_particles)=0;
  
  /*
  virtual void weight_for_adapting_sequence(const Index* index,
                                            Particles &current_particles,
                                            const Parameters &conditioned_on_parameters)=0;
  */
  
  virtual MoveOutput* subsample_move(RandomNumberGenerator &rng,
                                     Particle &particle)=0;
  
  /*
  virtual MoveOutput* subsample_move(RandomNumberGenerator &rng,
                                     Particle &particle,
                                     const Parameters &conditioned_on_parameters)=0;
  */
  
  //virtual void specific_subsample_move(RandomNumberGenerator &rng,
  //                                     Particle* new_particle,
  //                                     Particle &particle,
  //                                     const Parameters &conditioned_on_parameters)=0;
  
  virtual void subsample_weight_for_adapting_sequence(const Index* index,
                                                      Particles &current_particles)=0;
  
  /*
  virtual void subsample_weight_for_adapting_sequence(const Index* index,
                                                      Particles &current_particles,
                                                      const Parameters &conditioned_on_parameters)=0;
  */
  
  //EvaluateLogDistributionPtr evaluate_log_proposal;
  
protected:
  
  friend SequentialSMCWorker;
  friend class SimulateWorker;
  friend class ConditionalSimulateWorker;
  friend class SubsampleSimulateWorker;
  friend class SubsampleConditionalSimulateWorker;
  friend class MoveWorker;
  friend class SubsampleMoveWorker;
  
  bool proposal_is_evaluated;
  
  // stored here
  ParticleSimulator* particle_simulator;
  
  virtual SMCOutput* specific_run()=0;
  //virtual SMCOutput* specific_run(const std::string &directory_name)=0;

  void simulate_proposal(SMCOutput* simulation);
                                        
  virtual SMCOutput* specific_run(const Parameters &conditioned_on_parameters)=0;
  //virtual SMCOutput* specific_run(const std::string &directory_name,
  //                                const Parameters &conditioned_on_parameters)=0;
  
  void simulate_proposal(SMCOutput* simulation,
                         const Parameters &conditioned_on_parameters);
  
  // Do what needs to be done for each particlar type of SMC (might iterate, might not).
  SMCOutput* initialise_smc();
  virtual SMCOutput* specific_initialise_smc()=0;
  virtual void simulate_smc(SMCOutput* simulation)=0;
  //virtual void evaluate_smc(SMCOutput* simulation,const Parameters &conditioned_on_parameters)=0;
  virtual void evaluate_smc(SMCOutput* simulation)=0;
  virtual void evaluate_smcfixed_part_smc(SMCOutput* simulation)=0;
  virtual void evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation)=0;
  
  SMCOutput* initialise_smc(const Parameters &parameters);
  virtual SMCOutput* specific_initialise_smc(const Parameters &parameters)=0;
  virtual void simulate_smc(SMCOutput* simulation,
                            const Parameters &parameters)=0;
  //virtual void evaluate_smc(SMCOutput* simulation,const Parameters &conditioned_on_parameters)=0;
  virtual void evaluate_smc(SMCOutput* simulation,
                            const Parameters &parameters)=0;
  virtual void evaluate_smcfixed_part_smc(SMCOutput* simulation,
                                          const Parameters &parameters)=0;
  virtual void evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation,
                                                            const Parameters &parameters)=0;
  
  virtual void subsample_simulate_smc(SMCOutput* simulation)=0;
  virtual void subsample_evaluate_smc(SMCOutput* simulation)=0;
  virtual void subsample_evaluate_smcfixed_part_smc(SMCOutput* simulation)=0;
  virtual void subsample_evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation)=0;
  
  virtual void subsample_simulate_smc(SMCOutput* simulation,
                                      const Parameters &parameters)=0;
  virtual void subsample_evaluate_smc(SMCOutput* simulation,
                                      const Parameters &parameters)=0;
  virtual void subsample_evaluate_smcfixed_part_smc(SMCOutput* simulation,
                                                    const Parameters &parameters)=0;
  virtual void subsample_evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation,
                                                                      const Parameters &parameters)=0;
  
  void setup_default_ancestor_variables();
  
  //virtual void smc_simulate(SMCOutput* current_state)=0;
  //virtual void smc_weight(SMCOutput* current_state)=0;

  void make_copy(const SMC &another);

  //Particles is_step() const;

  friend SMCWorker;
  friend SMCOutput;
  // Stored here.
  SMCWorker* the_worker;
  // Stored here.
  SMCCriterion* resampling_criterion;
  
  // Stored here.
  //Transform* transform;
  //bool store_raw;
  //bool store_transformed;
  
  Sequencer sequencer;
  // not stored here
  //Parameters* sequencer_parameters;
  bool sequencer_limit_is_fixed;
  size_t number_of_particles;
  size_t lag;
  size_t lag_proposed;
  
  //SummaryStatisticsPtr summary_statistics;
  
  bool proposed_particles_inputted;
  std::vector<Parameters> initial_particles; // not needed, unless initial values provided
  arma::colvec log_probabilities_of_initial_values;
  
  std::vector<std::string> vector_variables;
  std::vector<std::string> any_variables;
  
  std::vector<size_t> vector_variable_sizes;
  
  std::vector<size_t> default_ancestor_variables;
  
  std::string results_name;
  
  bool initialised;
  
  std::ofstream log_likelihood_file_stream;
  std::ofstream time_file_stream;
  std::ofstream vector_variables_file_stream;
  std::ofstream vector_variable_sizes_file_stream;
  std::ofstream incremental_log_likelihood_file_stream;
  std::ofstream resampled_file_stream;
  std::ofstream ess_file_stream;
  std::ofstream schedule_parameters_file_stream;
  std::ofstream vector_points_file_stream;
  std::ofstream any_points_file_stream; // should be one for each member of Parameters
  std::ofstream normalised_weights_file_stream;
  std::ofstream unnormalised_weights_file_stream;
  
  //Parameters single_particle_is_step() const;

  //virtual void smc_step()=0;

  //virtual void weight_update()=0;

  //virtual double single_particle_weight_update() const=0;

  // Stored here.
  //SMCOutput* output;

};

#endif
