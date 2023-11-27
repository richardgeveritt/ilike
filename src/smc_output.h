#ifndef SMCOUTPUT_H
#define SMCOUTPUT_H

//#include <Rcpp.h>

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <deque>
#include <chrono>

#include "likelihood_estimator_output.h"
#include "particles.h"
#include "sequencer.h"

class SMC;
class ImportanceSampler;
class SMCMCMCMove;
class SMCMarginal;
class SMCGeneric;
class ParticleFilter;

class SMCOutput : public LikelihoodEstimatorOutput
{

public:

  SMCOutput();
  virtual ~SMCOutput();

  SMCOutput(SMC* estimator_in,
            size_t lag_in,
            size_t lag_proposed_in,
            const std::string &results_name_in);

  SMCOutput(const SMCOutput &another);
  void operator=(const SMCOutput &another);
  LikelihoodEstimatorOutput* duplicate() const;
  SMCOutput* smc_duplicate() const;

  // Additional functions compared to some other LikelihoodEstimators, since we often run SMC without it simply being to estimate a likelihood.
  void simulate();
  void evaluate_smcfixed_part();
  void evaluate_smcadaptive_part_given_smcfixed();
  
  void simulate(const Parameters &parameters);
  void evaluate_smcfixed_part(const Parameters &parameters);
  void evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters);
  
  void subsample_simulate();
  void subsample_evaluate_smcfixed_part();
  void subsample_evaluate_smcadaptive_part_given_smcfixed();
  
  void subsample_simulate(const Parameters &parameters);
  void subsample_evaluate_smcfixed_part(const Parameters &parameters);
  void subsample_evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters);
  
  Particles back() const;
  Particles& back();
  
  std::deque<Particles>::iterator end();
  std::deque<Particles>::const_iterator end() const;
  
  double calculate_latest_log_normalising_constant_ratio();
  
  arma::mat get_gradient_of_log(const std::string &variable,
                                const Parameters &x);
  arma::mat subsample_get_gradient_of_log(const std::string &variable,
                                          const Parameters &x);
  
  Particles* add_particles();
  Particles* add_particles(Particles* most_recent_particles);
  void add_proposed_particles(const Particles &latest_proposals);
  //void initialise_unnormalised_log_incremental_weights(const arma::colvec &latest_unnormalised_log_incremental_weights);
  //void initialise_next_step();
  void update_weights(const arma::colvec &latest_unnormalised_log_incremental_weights);
  void normalise_and_resample_weights();
  void resample();
  void mcmc_move();
  
  LikelihoodEstimator* get_likelihood_estimator() const;
  
  size_t number_of_smc_iterations() const;
  
  void increment_smc_iteration();
  void decrement_smc_iteration();
  
  void set_time();

  void print(std::ostream &os) const;
  
  double log_likelihood_pre_last_step;
  
  void close_ofstreams();

  void forget_you_were_already_written_to_file();
  
  std::deque<double> times;
  
  std::deque<double> llhds;

protected:
  
  //friend ImportanceSampler;
  //friend SMCMCMCMove;
  //friend SMC;
  //friend SMCMarginal;
  //friend SMCGeneric;
  //friend Sequencer;
  //friend ParticleFilter;
  
  // Stored in ModelAndAlgorithm or in main.
  SMC* estimator;

  std::deque<Particles> all_particles;

  std::deque<Particles> all_proposed;

  size_t lag;

  size_t lag_proposed;
  
  size_t smc_iteration;
  
  int iteration_written_to_file;
  
  std::string results_name;
  
  std::chrono::high_resolution_clock::time_point start_time;
  
  void write_to_file(const std::string &directory_name,
                     const std::string &index = "");
  
  void close_ofstreams(size_t deque_index);

  void make_copy(const SMCOutput &another);

};

#endif
