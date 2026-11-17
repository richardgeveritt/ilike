#ifndef SMCOUTPUT_H
#define SMCOUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <deque>
#include <chrono>

#include "likelihood_estimator_output.h"
#include "particles.h"
#include "sequencer.h"

namespace ilike
{

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
  int adaptive_summary_size() const;
  void adaptive_summary(double *values) const;
  double smcadaptive_from_summary(const double *values,
                                  const Parameters &parameters) const;
  bool repoint_estimator(LikelihoodEstimator* estimator_in);
  
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
  void update_weights(const arma::colvec &latest_unnormalised_log_incremental_weights);
  void normalise_and_resample_weights();
  void resample();
  void mcmc_move();
  
  LikelihoodEstimator* get_likelihood_estimator() const;
  
  size_t number_of_smc_iterations() const;
  
  void increment_smc_iteration();
  void decrement_smc_iteration();
  
  void set_time();
  void set_llhd(double llhd_in);
  
  void set_time_and_reset_start();
  
  void print(std::ostream &os) const;
  
  double log_likelihood_pre_last_step;
  
  void close_ofstreams();
  
  void forget_you_were_already_written_to_file();
  
  void terminate();
  
  std::deque<double> times;
  
  std::deque<double> llhds;
  
protected:

  // Stored in ModelAndAlgorithm or in main.
  SMC* estimator;
  
  std::deque<Particles> all_particles;

  // The stored generations of a nested run in compact form (see
  // SMC::compact_stored_output): once compacted, all_particles is empty and
  // these are what is written.
  struct CompactGeneration
  {
    size_t iteration;
    arma::mat points_transposed;   // one column per point
    arma::colvec unnormalised_log_weights;
    arma::colvec normalised_log_weights;
    std::vector<size_t> ancestors; // empty if the previous generation was not stored
    std::vector<double> output_lengths;
    double incremental_log_likelihood;
    double ess;
    int resampled;
    std::string schedule_parameters;
  };
  std::vector<CompactGeneration> compact_generations;
  bool compacted = false;
  void compact();
  void write_compact_nested(const std::string &dir_name);
  
  std::deque<Particles> all_proposed;
  
  size_t lag;
  size_t output_lag;
  
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
}

#endif
