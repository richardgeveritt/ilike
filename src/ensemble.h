#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include "particle.h"

class SequentialEnsembleKalmanWorker;
class PackingInstructions;
class EnsembleFactors;
class AdjustmentEnsembleShifter;
class MoveOutput;
class GaussianSMCAdaptor;
class EnsembleKalmanMFDS;

class Ensemble
{
public:

	Ensemble();
	virtual ~Ensemble();

	// Everything you need to copy the class.
	Ensemble(const Ensemble &another);
	Ensemble* duplicate() const;
	void make_copy(const Ensemble &another);
	void operator=(const Ensemble &another);
  
  void reserve(size_t number_of_ensemble_members_in);
  //void push_back(const Parameters &parameters_in);
  void push_back(const Particle &ensemble_member_in);
  void push_back(MoveOutput* move_output_in);
  
  MoveOutput* operator[](const size_t &i);
  MoveOutput* operator[](const size_t &i) const;
  
  MoveOutput* back();
  MoveOutput* back() const;
  
  size_t size() const;
  
  void find_measurement_covariances();
  void find_measurement_covariances(const Parameters &conditioned_on_parameters);
  
  void set_temperature(double temperature_in);
  
  void set_previous_target_evaluated_to_target_evaluated();
  void subsample_set_previous_target_evaluated_to_target_evaluated();

protected: // Things that can be accessed in this class and subclasses.
  
  friend SequentialEnsembleKalmanWorker;
  friend AdjustmentEnsembleShifter;
  friend GaussianSMCAdaptor;
  friend EnsembleKalmanMFDS;
  
  std::vector<MoveOutput*> members;
  std::vector<arma::colvec> partially_packed_members_col;
  std::vector<arma::rowvec> partially_packed_members_row;
  arma::mat packed_members;
  
  // inner vector is different measurements, outer is realisations from different members of the ensemble
  std::vector< std::vector<arma::rowvec> > partially_packed_measurement_states;
  //std::vector< std::vector<arma::rowvec> > partially_packed_measurement_random_shifts;
  
  std::vector<arma::mat> packed_measurement_states;
  
  // vectors over different measurements
  // first two of these could be stored in EnsembleFactors - a bunch of the work that relies on storing vectors of these things (where we don't know that we really have a vector of factors) could be done in that class.
  std::vector<arma::mat> Cxys;
  std::vector<arma::mat> Cyys; // need to rename since this is the raw estimate, but it is not exactly Cyy (needs +\Lambda)
  arma::mat inv_Cxx; // not calculated unless we can help it (have intractable llhds)
  
  // vector of different measurements
  // should be vector, but not written yet
  //arma::colvec measurements;
  
  // not stored here
  PackingInstructions* packing_instructions;
  EnsembleFactors* ensemble_factors;

private: // Things that can be accessed only by this class.

};

#endif
