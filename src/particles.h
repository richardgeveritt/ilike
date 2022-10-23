#ifndef PARTICLES_H
#define PARTICLES_H

//#include <RcppParallel.h>
#include <RcppArmadillo.h>
#include <iostream>

using namespace Rcpp;
//using namespace RcppParallel;

#include "particle.h"
#include "distributions.h"

class SMC;
class SequentialSMCWorker;
class MoveOutput;

class Particles
{
public:

  Particles(void);
  Particles(size_t number_of_particles_in);
  Particles(const std::vector< MoveOutput* > &particles_in);
  Particles(const std::vector<Parameters> &initial_values_in,
            const arma::colvec &log_probabilities_of_initial_values_in);

  virtual ~Particles(void);

  Particles(const Particles &another);
  void operator=(const Particles &another);
  
  void reserve(size_t number_of_particles_in);
  //void push_back(const Parameters &parameters_in);
  void push_back(const Particle &particle_in);
  void push_back(MoveOutput* move_output_in);
  
  void simulate_resampling_variables(RandomNumberGenerator &rng);

  //arma::colvec get_resampling_variables() const;
  //void set_ancestor_variables(const std::vector<size_t> &ancestor_variables_in);
  
  void initialise_weights();
  void update_weights(const arma::colvec &latest_unnormalised_log_incremental_weights);
  void normalise_weights();
  void resample();
  
  void set_previous_target_evaluated_to_target_evaluated();
  void subsample_set_previous_target_evaluated_to_target_evaluated();
  
  arma::mat get_most_recent_matrix_particles() const;
  
  //virtual List& operator[](int index)=0;
  //virtual List operator[](int index) const=0;

  friend std::ostream& operator<<(std::ostream& os, const Particles &p);
  
  // Access only elements in std::vector<Particle>.
  //double& operator[](const size_t &i);
  //double operator[](const size_t &i) const;
  
  MoveOutput* operator[](const size_t &i);
  MoveOutput* operator[](const size_t &i) const;
  
  MoveOutput* back();
  MoveOutput* back() const;
  
  size_t size() const;
  
  // Outer vector is over particles, inner is stored here
  std::vector< MoveOutput* > particles;
  
  arma::colvec unnormalised_log_weights;
  
  arma::colvec incremental_log_weights;

  arma::colvec normalised_log_weights;
  
  arma::colvec previous_normalised_log_weights;

  double log_normalising_constant_ratio;
  
  // Also store anything that needs to be stored in a different structure across all particles.
  arma::colvec resampling_variables;
  std::vector<size_t> ancestor_variables;

protected:
  
  //friend SequentialSMCWorker;

  void make_copy(const Particles &another);

  // Stored here. // moved to model_and_algorithm
  //std::vector<LikelihoodEstimator*> likelihood_estimators;

};

#endif
