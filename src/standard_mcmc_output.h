#ifndef STANDARDMCMCOUTPUT_H
#define STANDARDMCMCOUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "move_output.h"
#include "parameters.h"
#include "ilike_header.h"

class MCMCTermination;
class MCMC;

class StandardMCMCOutput : public MoveOutput
{

public:

  //StandardMCMCOutput(const Parameters &parameter_in);
  StandardMCMCOutput();
  
  StandardMCMCOutput(MCMCTermination* termination_in);

  virtual ~StandardMCMCOutput();

  StandardMCMCOutput(const StandardMCMCOutput &another);

  void operator=(const StandardMCMCOutput &another);
  //MoveOutput* duplicate() const;
  
  void push_back(const Particle &particle_in);

  Particle& back();
  Particle back() const;
  
  std::vector<Parameters> get_vector_of_parameters() const;
  
  void write_vector_points(const std::vector<std::string> &variables,
                           std::ofstream &file_stream,
                           std::shared_ptr<Transform> transform) const;
  void write_any_points(const std::vector<std::string> &variables,
                        std::ofstream &file_stream) const;
  
  void write_factors(const std::string &directory_name,
                     const std::string &index) const;
  void write_ensemble_factors(const std::string &directory_name,
                              const std::string &index) const;
  
  size_t length() const;
  
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
  
protected:
  
  //virtual void specific_mcmc_adapt()=0;
  
  size_t iteration_counter;
  
  // Stored here.
  MCMCTermination* termination;
  
  void make_copy(const StandardMCMCOutput &another);
  
  std::deque<Particle> output;

};

#endif
