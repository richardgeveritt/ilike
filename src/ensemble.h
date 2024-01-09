#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include "particle.h"

class SequentialEnsembleKalmanWorker;
class PackingInstructions;
class EnsembleFactors;
class AdjustmentEnsembleShifter;
class SquareRootEnsembleShifter;
class MoveOutput;
class GaussianSMCAdaptor;
class EnsembleKalmanMFDS;
class VectorEnsembleFactors;
class HMMEnsembleFactors;
class ESSSMCCriterion;

class Ensemble
{
public:

	Ensemble();
	virtual ~Ensemble();
  
  Ensemble(EnsembleFactors* ensemble_factors_in);
  Ensemble(std::vector<Parameters> &initial_values_in,
           EnsembleFactors* factors_in,
           const Parameters &conditioned_on_parameters);
  
  Ensemble(std::vector<Parameters> &initial_values_in,
           EnsembleFactors* factors_in,
           const Parameters &conditioned_on_parameters,
           const Parameters &sequencer_parameters);

	// Everything you need to copy the class.
	Ensemble(const Ensemble &another);
  Ensemble& operator=(const Ensemble &another);
  
  Ensemble(Ensemble &&another);
  Ensemble& operator=(Ensemble &&another);
  
	Ensemble* duplicate() const;
	
  void setup(std::vector<Parameters> &initial_values_in,
             EnsembleFactors* factors_in,
             const Parameters &conditioned_on_parameters);
  
  void setup(std::vector<Parameters> &initial_values_in,
             EnsembleFactors* factors_in,
             const Parameters &conditioned_on_parameters,
             const Parameters &sequencer_parameters);
  
  void reserve(size_t number_of_ensemble_members_in);
  //void push_back(const Parameters &parameters_in);
  void push_back(Particle &&ensemble_member_in);
  void push_back(MoveOutput* move_output_in);
  void push_back(Parameters &&parameters_in,
                 EnsembleFactors* factors_in);
  
  Particle* add_ensemble_member();
  
  MoveOutput* operator[](const size_t &i);
  MoveOutput* operator[](const size_t &i) const;
  
  MoveOutput* back();
  MoveOutput* back() const;
  
  size_t size() const;
  
  arma::rowvec get_output_lengths() const;
  
  void find_measurement_covariances();
  void find_measurement_covariances(const Parameters &conditioned_on_parameters);
  
  void set_temperature(double temperature_in);
  double get_inverse_incremental_temperature() const;
  void update_weights(const arma::colvec &latest_unnormalised_log_incremental_weights);
  
  double calculate_log_normalising_constant();
  double calculate_inversion_log_normalising_constant(double inverse_incremental_temperature);
  
  void set_previous_target_evaluated_to_target_evaluated();
  void subsample_set_previous_target_evaluated_to_target_evaluated();
  
  void precompute_gaussian_covariance(double inverse_incremental_temperature);
  
  arma::mat get_packed_members() const;
  
  double log_normalising_constant_ratio;
  
  // stored here
  std::vector<MoveOutput*> members;
  std::vector<MoveOutput*> predicted_members;
  
  // only here for output
  Parameters schedule_parameters;
  
  double ess;
  
  void close_ofstreams();

protected: // Things that can be accessed in this class and subclasses.
  
  void make_copy(const Ensemble &another);
  void make_copy(Ensemble &&another);
  
  friend SequentialEnsembleKalmanWorker;
  friend AdjustmentEnsembleShifter;
  friend SquareRootEnsembleShifter;
  //friend GaussianSMCAdaptor;
  friend EnsembleKalmanMFDS;
  friend VectorEnsembleFactors;
  friend HMMEnsembleFactors;
  friend ESSSMCCriterion;
  
  std::vector<arma::colvec> partially_packed_members_col;
  std::vector<arma::colvec> partially_packed_predicted_members_col;
  std::vector<arma::rowvec> partially_packed_members_row;
  arma::mat packed_members;
  
  // inner vector is different measurements, outer is realisations from different members of the ensemble
  std::vector< std::vector<arma::rowvec> > partially_packed_measurement_states;
  //std::vector< std::vector<arma::rowvec> > partially_packed_measurement_random_shifts;
  
  std::vector<arma::mat> packed_measurement_states;
  
  // vectors over different measurements
  // first two of these could be stored in EnsembleFactors - a bunch of the work that relies on storing vectors of these things (where we don't know that we really have a vector of factors) could be done in that class.
  std::vector<arma::mat> Cxys;
  std::vector<arma::colvec> myys;
  std::vector<arma::mat> Cyys; // need to rename since this is the raw estimate, but it is not exactly Cyy (needs +\Lambda)
  std::vector<arma::mat> kalman_gains;
  arma::mat inv_Cxx; // not calculated unless we can help it (have intractable llhds)
  
  // vector of different measurements
  // should be vector, but not written yet
  //arma::colvec measurements;
  
  arma::colvec unnormalised_log_weights;
  
  // not stored here
  //PackingInstructions* packing_instructions;
  EnsembleFactors* ensemble_factors;

private: // Things that can be accessed only by this class.

};

#endif
