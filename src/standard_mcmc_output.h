#ifndef STANDARDMCMCOUTPUT_H
#define STANDARDMCMCOUTPUT_H

#include <RcppArmadillo.h>
#include <memory>
using namespace Rcpp;

#include <vector>
#include "move_output.h"
#include "particle.h"
#include "sample_store.h"
#include "parameters.h"
#include "ilike_header.h"

namespace ilike
{

class MCMCTermination;
class MCMC;

class StandardMCMCOutput : public MoveOutput
{
  
public:
  
  StandardMCMCOutput();
  
  StandardMCMCOutput(MCMCTermination* termination_in);
  
  virtual ~StandardMCMCOutput();
  
  StandardMCMCOutput(const StandardMCMCOutput &another);
  
  void operator=(const StandardMCMCOutput &another);
  
  void push_back(const Particle &particle_in);

  // Store an iteration (subject to thin) without changing the current state.
  void record(const Particle &particle_in);

  // Set the current state, taking it over.
  void set_current(Particle &&particle_in);

  void set_thin(size_t thin_in);

  void mcmc_adapt(const Particle &current_particle);
  
  Particle& back();
  Particle back() const;
  
  std::vector<Parameters> get_vector_of_parameters() const;
  
  void write_vector_points(const std::vector<std::string> &variables,
                           std::ofstream &file_stream,
                           std::shared_ptr<Transform> transform) const;
  void write_any_points(const std::vector<std::string> &variables,
                        std::ofstream &file_stream) const;

  size_t fill_vector_points_transposed(const std::vector<std::string> &variables,
                                       arma::mat &destination,
                                       size_t first_column) const;

  arma::mat get_matrix_of_vector_points(const std::vector<std::string> &variables,
                                        std::shared_ptr<Transform> transform) const;
  
  void write_factors(const std::string &directory_name,
                     const std::string &index) const;
  void write_ensemble_factors(const std::string &directory_name,
                              const std::string &index) const;
  
  size_t length() const;
  bool repoint_to(const Factors* factors);
  
  void close_ofstreams();
  
  size_t* get_iteration_counter_pointer();
  
  void increment_counter();
  
  void reset_counter();
  
  void mcmc_adapt();
  
  Particle move(RandomNumberGenerator &rng,
                const Particle &particle) const;
  
  Particle subsample_move(RandomNumberGenerator &rng,
                          const Particle &particle) const;
  
  bool terminate() const;
  
  virtual MCMC* get_mcmc()=0;
  virtual const MCMC* get_mcmc() const=0;
  
  Parameters get_current_algorithm_parameters() const;

  // Adaptation output (dev/plan.md 3.2): from start_recording_adaptation on,
  // each stored iteration also stores the move's adaptation trace;
  // finish_adaptation stores its final state.
  void start_recording_adaptation();
  void finish_adaptation();
  size_t adaptation_trace_columns() const;
  size_t adaptation_final_columns() const;
  size_t fill_adaptation_trace_transposed(arma::mat &destination,
                                          size_t first_column) const;
  void fill_adaptation_final(arma::mat &destination,
                             size_t column) const;
  void adaptation_names(std::vector<std::string> &trace_names,
                        std::vector<std::string> &final_names) const;
  
protected:

  bool records_adaptation = false;
  size_t adaptation_columns = 0;
  std::vector<double> adaptation_trace; // stored iterations, row after row
  std::vector<double> adaptation_final;

  size_t iteration_counter;
  
  // Stored here.
  // owned (dev/plan.md 3.5, CH3)
  std::unique_ptr<MCMCTermination> termination;
  
  void make_copy(const StandardMCMCOutput &another);
  
  // The current state of the chain: what back() returns.
  Particle current;

  // Stored iterations (dev/plan.md 2.5, "Stage 2 - concrete design", item
  // 2): their parameter values and, as extra columns, one compact record per
  // factor (LikelihoodEstimatorOutput::compact_value). A stored MCMC
  // iteration then costs 8 (d + number of factors) bytes instead of a
  // Particle (about 1.7 KB at d = 1).
  SampleStore store;

  // Full copies of the stored iterations, used instead of the store when an
  // iteration cannot be stored compactly: ensemble factor variables, or a
  // factor output without a compact record (density/synthetic likelihood,
  // annealed, nested algorithms, ...). Decided at the first stored
  // iteration (rich_storage).
  std::deque<Particle> output;
  bool storage_decided;
  bool rich_storage;

  // Store every thin-th iteration only (the current state is always kept).
  size_t thin;
  size_t number_of_iterations_seen;

  std::vector<double> compact_scratch;
  
  std::deque<Parameters> algorithm_parameters;
  
};
}

#endif
