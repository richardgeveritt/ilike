#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include "particle.h"

namespace ilike
{
  /**
   * @file ensemble.h
   * @brief Defines the SequentialEnsembleKalmanWorker class.
   *
   * An Ensemble Kalman sequential worker implementation. Performs ensemble-based data assimilation using Kalman-type updates.
   *
   * @namespace ilike
   * @class SequentialEnsembleKalmanWorker
   * @brief The sequential ensemble kalman worker class.
   */


class SequentialEnsembleKalmanWorker;
class PackingInstructions;
class EnsembleFactors;
class AdjustmentEnsembleShifter;
class SquareRootEnsembleShifter;
class MoveOutput;
class GaussianSMCAdaptor;
class EnsembleKalmanMFDS;
class EnsembleKalmanInversion;
class VectorEnsembleFactors;
class HMMEnsembleFactors;
class ESSSMCCriterion;

class Ensemble
{
public:
  
  /**
   * @brief Performs the ensemble operation.
   */
  Ensemble();
  /**
   * @brief Performs the ~ensemble operation.
   */
  virtual ~Ensemble();
  
  /**
   * @brief Performs the ensemble operation.
   *
   * @param ensemble_factors_in The ensemble factors.
   */
  Ensemble(EnsembleFactors* ensemble_factors_in);
  Ensemble(std::vector<Parameters> &initial_values_in,
           EnsembleFactors* factors_in,
           const Parameters &conditioned_on_parameters);
  
  Ensemble(std::vector<Parameters> &initial_values_in,
           EnsembleFactors* factors_in,
           const Parameters &conditioned_on_parameters,
           const Parameters &sequencer_parameters);
  
  // Everything you need to copy the class.
  /**
   * @brief Performs the ensemble operation.
   *
   * @param another The SequentialEnsembleKalmanWorker instance to copy from.
   */
  Ensemble(const Ensemble &another);
  /**
   * @brief Assignment operator for SequentialEnsembleKalmanWorker.
   *
   * @param another The SequentialEnsembleKalmanWorker instance to copy from.
   */
  Ensemble& operator=(const Ensemble &another);
  
  /**
   * @brief Performs the ensemble operation.
   *
   * @param another The SequentialEnsembleKalmanWorker instance to copy from.
   */
  Ensemble(Ensemble &&another);
  /**
   * @brief Assignment operator for SequentialEnsembleKalmanWorker.
   *
   * @param another The SequentialEnsembleKalmanWorker instance to copy from.
   */
  Ensemble& operator=(Ensemble &&another);
  
  /**
   * @brief Creates a deep copy of this SequentialEnsembleKalmanWorker object.
   *
   * @return The result.
   */
  Ensemble* duplicate() const;
  
  void setup(std::vector<Parameters> &initial_values_in,
             EnsembleFactors* factors_in,
             const Parameters &conditioned_on_parameters);
  
  void setup(std::vector<Parameters> &initial_values_in,
             EnsembleFactors* factors_in,
             const Parameters &conditioned_on_parameters,
             const Parameters &sequencer_parameters);
  
  /**
   * @brief Performs the reserve operation.
   *
   * @param number_of_ensemble_members_in The number of ensemble members.
   */
  void reserve(size_t number_of_ensemble_members_in);
  //void push_back(const Parameters &parameters_in);
  /**
   * @brief Performs the push back operation.
   *
   * @param ensemble_member_in The ensemble member.
   */
  void push_back(Particle &&ensemble_member_in);
  /**
   * @brief Performs the push back operation.
   *
   * @param move_output_in The move output.
   */
  void push_back(MoveOutput* move_output_in);
  void push_back(Parameters &&parameters_in,
                 EnsembleFactors* factors_in);
  
  /**
   * @brief Performs the add ensemble member operation.
   *
   * @return The result.
   */
  Particle* add_ensemble_member();
  
  /**
   * @brief Performs the operator[] operation.
   *
   * @param i The i.
   *
   * @return The result.
   */
  MoveOutput* operator[](const size_t &i);
  /**
   * @brief Performs the operator[] operation.
   *
   * @param i The i.
   *
   * @return The result.
   */
  MoveOutput* operator[](const size_t &i) const;
  
  /**
   * @brief Performs the back operation.
   *
   * @return The result.
   */
  MoveOutput* back();
  /**
   * @brief Performs the back operation.
   *
   * @return The result.
   */
  MoveOutput* back() const;
  
  /**
   * @brief Performs the size operation.
   *
   * @return The result.
   */
  size_t size() const;
  
  /**
   * @brief Returns the output lengths.
   *
   * @return The result.
   */
  arma::rowvec get_output_lengths() const;
  
  /**
   * @brief Finds measurement covariances.
   */
  void find_measurement_covariances();
  /**
   * @brief Finds measurement covariances.
   *
   * @param conditioned_on_parameters The conditioned on parameters.
   */
  void find_measurement_covariances(const Parameters &conditioned_on_parameters);
  
  /**
   * @brief Sets the temperature.
   *
   * @param temperature_in The temperature.
   */
  void set_temperature(double temperature_in);
  /**
   * @brief Returns the inverse incremental temperature.
   *
   * @return The result.
   */
  double get_inverse_incremental_temperature() const;
  /**
   * @brief Updates weights.
   *
   * @param latest_unnormalised_log_incremental_weights The latest unnormalised log incremental weights.
   */
  void update_weights(const arma::colvec &latest_unnormalised_log_incremental_weights);
  
  /**
   * @brief Performs the calculate log normalising constant operation.
   *
   * @return The result.
   */
  double calculate_log_normalising_constant();
  /**
   * @brief Performs the calculate mc inversion log normalising constant operation.
   *
   * @param unnormalised_log_weights The unnormalised log weights.
   * @param inverse_incremental_temperature The inverse incremental temperature.
   *
   * @return The result.
   */
  double calculate_mc_inversion_log_normalising_constant(const arma::colvec &unnormalised_log_weights, double inverse_incremental_temperature);
  /**
   * @brief Performs the calculate inversion log normalising constant operation.
   *
   * @param inverse_incremental_temperature The inverse incremental temperature.
   *
   * @return The result.
   */
  double calculate_inversion_log_normalising_constant(double inverse_incremental_temperature);
  /**
   * @brief Performs the calculate unbiased inversion log normalising constant operation.
   *
   * @param inverse_incremental_temperature The inverse incremental temperature.
   *
   * @return The result.
   */
  double calculate_unbiased_inversion_log_normalising_constant(double inverse_incremental_temperature);
  void calculate_path1_inversion_log_normalising_constant(std::vector<double> &log_measurement_likelihood_means,
                                                          double temperature,
                                                          double multiplier);
  void calculate_path2_inversion_log_normalising_constant(std::vector<double> &log_measurement_likelihood_means,
                                                          std::vector<double> &log_measurement_likelihood_variances);
  
  /**
   * @brief Performs the calculate kalman gains operation.
   *
   * @param inverse_incremental_temperature The inverse incremental temperature.
   */
  void calculate_kalman_gains(double inverse_incremental_temperature);
  
  /**
   * @brief Sets the previous target evaluated to target evaluated.
   */
  void set_previous_target_evaluated_to_target_evaluated();
  /**
   * @brief Performs the subsample set previous target evaluated to target evaluated operation.
   */
  void subsample_set_previous_target_evaluated_to_target_evaluated();
  
  /**
   * @brief Performs the precompute gaussian covariance operation.
   *
   * @param inverse_incremental_temperature The inverse incremental temperature.
   */
  void precompute_gaussian_covariance(double inverse_incremental_temperature);
  
  /**
   * @brief Returns the packed members.
   *
   * @return The result.
   */
  arma::mat get_packed_members() const;
  
  double log_normalising_constant_ratio;
  
  // stored here
  std::vector<MoveOutput*> members;
  std::vector<MoveOutput*> predicted_members;
  
  // only here for output
  Parameters schedule_parameters;
  
  double ess;
  
  // used in path sampling, otherwise empty
  std::vector<double> log_measurement_likelihood_means;
  std::vector<double> log_measurement_likelihood_variances;
  
  /**
   * @brief Closes any open file streams.
   */
  void close_ofstreams();
  
protected: // Things that can be accessed in this class and subclasses.
  
  /**
   * @brief Copies the state of another SequentialEnsembleKalmanWorker into this object.
   *
   * @param another The SequentialEnsembleKalmanWorker instance to copy from.
   */
  void make_copy(const Ensemble &another);
  /**
   * @brief Copies the state of another SequentialEnsembleKalmanWorker into this object.
   *
   * @param another The SequentialEnsembleKalmanWorker instance to copy from.
   */
  void make_copy(Ensemble &&another);
  
  friend SequentialEnsembleKalmanWorker;
  friend AdjustmentEnsembleShifter;
  friend SquareRootEnsembleShifter;
  //friend GaussianSMCAdaptor;
  friend EnsembleKalmanMFDS;
  friend EnsembleKalmanInversion;
  friend VectorEnsembleFactors;
  friend HMMEnsembleFactors;
  friend ESSSMCCriterion;
  
  /** @brief The partially packed members col. */
  std::vector<arma::colvec> partially_packed_members_col;
  /** @brief The partially packed predicted members col. */
  std::vector<arma::colvec> partially_packed_predicted_members_col;
  /** @brief The partially packed members row. */
  std::vector<arma::rowvec> partially_packed_members_row;
  /** @brief The packed members. */
  arma::mat packed_members;
  
  // inner vector is different measurements, outer is realisations from different members of the ensemble
  /** @brief The partially packed measurement states. */
  std::vector< std::vector<arma::rowvec> > partially_packed_measurement_states;
  //std::vector< std::vector<arma::rowvec> > partially_packed_measurement_random_shifts;
  
  /** @brief The packed measurement states. */
  std::vector<arma::mat> packed_measurement_states;
  
  // vectors over different measurements
  // first two of these could be stored in EnsembleFactors - a bunch of the work that relies on storing vectors of these things (where we don't know that we really have a vector of factors) could be done in that class.
  /** @brief The cxys. */
  std::vector<arma::mat> Cxys;
  /** @brief The myys. */
  std::vector<arma::colvec> myys;
  std::vector<arma::mat> Cyys; // need to rename since this is the raw estimate, but it is not exactly Cyy (needs +\Lambda)
  //std::vector<arma::mat> Cygivenxs;
  /** @brief The kalman gains. */
  std::vector<arma::mat> kalman_gains;
  
  /** @brief The inv sigma precomps. */
  std::vector<arma::mat> inv_sigma_precomps;
  /** @brief The log det precomps. */
  std::vector<double> log_det_precomps;
  
  arma::mat inv_Cxx; // not calculated unless we can help it (have intractable llhds)
  
  // vector of different measurements
  // should be vector, but not written yet
  //arma::colvec measurements;
  
  /** @brief The unnormalised log weights. */
  arma::colvec unnormalised_log_weights;
  
  // not stored here
  //PackingInstructions* packing_instructions;
  /** @brief The ensemble factors. */
  EnsembleFactors* ensemble_factors;
  
private: // Things that can be accessed only by this class.
  
};
}

#endif
