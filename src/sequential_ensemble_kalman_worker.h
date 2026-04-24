#ifndef SEQUENTIALENSEMBLEKALMANWORKER_H
#define SEQUENTIALENSEMBLEKALMANWORKER_H

#include "ensemble_kalman_worker.h"

#include <vector>

#include "ensemble_kalman.h"

namespace ilike
{
  /**
   * @file sequential_ensemble_kalman_worker.h
   * @brief Defines the SequentialEnsembleKalmanWorker class.
   *
   * An Ensemble Kalman sequential worker implementation. Performs ensemble-based data assimilation using Kalman-type updates.
   *
   * @namespace ilike
   * @class SequentialEnsembleKalmanWorker
   * @brief A sequential ensemble kalman worker derived from EnsembleKalmanWorker.
   */


class SequentialEnsembleKalmanWorker : public EnsembleKalmanWorker
{
public:
  
  /**
   * @brief Default constructor for SequentialEnsembleKalmanWorker.
   */
  SequentialEnsembleKalmanWorker();
  /**
   * @brief Destructor for SequentialEnsembleKalmanWorker.
   */
  virtual ~SequentialEnsembleKalmanWorker();
  
  /**
   * @brief Constructs a SequentialEnsembleKalmanWorker object.
   *
   * @param the_enk_in The the enk.
   */
  SequentialEnsembleKalmanWorker(EnsembleKalman* the_enk_in);
  
  /**
   * @brief Copy constructor for SequentialEnsembleKalmanWorker.
   *
   * @param another The SequentialEnsembleKalmanWorker instance to copy from.
   */
  SequentialEnsembleKalmanWorker(const SequentialEnsembleKalmanWorker &another);
  /**
   * @brief Assignment operator for SequentialEnsembleKalmanWorker.
   *
   * @param another The SequentialEnsembleKalmanWorker instance to copy from.
   */
  void operator=(const SequentialEnsembleKalmanWorker &another);
  /**
   * @brief Creates a deep copy of this SequentialEnsembleKalmanWorker object.
   *
   * @return The result.
   */
  EnsembleKalmanWorker* duplicate() const;
  
  void shift(Ensemble* ensemble,
             double inverse_incremental_temperature);
  
  /**
   * @brief Performs the pack operation.
   *
   * @param ensemble The ensemble.
   */
  void pack(Ensemble* ensemble);
  /**
   * @brief Performs the unpack operation.
   *
   * @param ensemble The ensemble.
   */
  void unpack(Ensemble* ensemble);
  /**
   * @brief Performs the unpack with predicted operation.
   *
   * @param ensemble The ensemble.
   */
  void unpack_with_predicted(Ensemble* ensemble);
  
  void weight(Ensemble* ensemble,
              const Index* index,
              double incremental_temperature);
  /*
   void weight(Ensemble* ensemble,
   const Index* index,
   double incremental_temperature,
   const Parameters &conditioned_on_parameters);
   */
  void subsample_weight(Ensemble* ensemble,
                        const Index* index,
                        double incremental_temperature);
  /*
   void subsample_weight(Ensemble* ensemble,
   const Index* index,
   double incremental_temperature,
   const Parameters &conditioned_on_parameters);
   */
  
  arma::colvec get_unnormalised_log_incremental_weights() const;
  
  /*
   void weight(EnsembleKalman &current_ensemble_kalman);
   void smcfixed_weight(EnsembleKalman &current_ensemble_kalman);
   void smcadaptive_given_smcfixed_weight(EnsembleKalman &current_ensemble_kalman);
   void smcadaptive_given_smcfixed_evaluate_target(EnsembleKalman &current_ensemble_kalman);
   void marginal_weight(EnsembleKalman &current_ensemble_kalman,
   EnsembleKalman &previous_ensemble_kalman,
   ProposalKernel* proposal_kernel);
   void generic_weight(EnsembleKalman &current_ensemble_kalman,
   EnsembleKalman &previous_ensemble_kalman,
   ProposalKernel* proposal_kernel,
   ProposalKernel* L_kernel);
   void pf_weight(EnsembleKalman &current_ensemble_kalman,
   EnsembleKalman &previous_ensemble_kalman,
   ProposalKernel* proposal_kernel);
   
   void weight(EnsembleKalman &current_ensemble_kalman,
   const Parameters &conditioned_on_parameters);
   void pf_initial_weight(EnsembleKalman &current_ensemble_kalman,
   const Parameters &conditioned_on_parameters);
   void smcfixed_weight(EnsembleKalman &current_ensemble_kalman,
   const Parameters &conditioned_on_parameters);
   void smcadaptive_given_smcfixed_weight(EnsembleKalman &current_ensemble_kalman,
   const Parameters &conditioned_on_parameters);
   void smcadaptive_given_smcfixed_evaluate_target(EnsembleKalman &current_ensemble_kalman,
   const Parameters &conditioned_on_parameters);
   void marginal_weight(EnsembleKalman &current_ensemble_kalman,
   EnsembleKalman &previous_ensemble_kalman,
   ProposalKernel* proposal_kernel,
   const Parameters &conditioned_on_parameters);
   void generic_weight(EnsembleKalman &current_ensemble_kalman,
   EnsembleKalman &previous_ensemble_kalman,
   ProposalKernel* proposal_kernel,
   ProposalKernel* L_kernel,
   const Parameters &conditioned_on_parameters);
   void pf_weight(EnsembleKalman &current_ensemble_kalman,
   EnsembleKalman &previous_ensemble_kalman,
   ProposalKernel* proposal_kernel,
   const Parameters &conditioned_on_parameters);
   
   void subsample_weight(EnsembleKalman &current_ensemble_kalman,
   const Parameters &conditioned_on_parameters);
   void subsample_smcfixed_weight(EnsembleKalman &current_ensemble_kalman,
   const Parameters &conditioned_on_parameters);
   void subsample_smcadaptive_given_smcfixed_weight(EnsembleKalman &current_ensemble_kalman,
   const Parameters &conditioned_on_parameters);
   void subsample_smcadaptive_given_smcfixed_evaluate_target(EnsembleKalman &current_ensemble_kalman,
   const Parameters &conditioned_on_parameters);
   void subsample_marginal_weight(EnsembleKalman &current_ensemble_kalman,
   EnsembleKalman &previous_ensemble_kalman,
   ProposalKernel* proposal_kernel,
   const Parameters &conditioned_on_parameters);
   void subsample_generic_weight(EnsembleKalman &current_ensemble_kalman,
   EnsembleKalman &previous_ensemble_kalman,
   ProposalKernel* proposal_kernel,
   ProposalKernel* L_kernel,
   const Parameters &conditioned_on_parameters);
   void subsample_pf_weight(EnsembleKalman &current_ensemble_kalman,
   EnsembleKalman &previous_ensemble_kalman,
   ProposalKernel* proposal_kernel,
   const Parameters &conditioned_on_parameters);
   
   //EnsembleKalman simulated_ensemble_kalman() const;
   //EnsembleKalman& simulated_ensemble_kalman();
   arma::colvec get_unnormalised_log_incremental_weights() const;
   */
  
protected:
  
  void specific_simulate(Ensemble* next_ensemble_kalman,
                         const Index* index);
  //void specific_simulate_and_weight(const Parameters &conditioned_on_parameters);
  
  void specific_simulate(Ensemble* next_ensemble_kalman,
                         const Index* index,
                         const Parameters &conditioned_on_parameters);
  
  void specific_move(Ensemble* next_particles,
                     Ensemble* current_particles);
  
  /*
   void specific_move(Ensemble* next_particles,
   Ensemble* current_particles,
   const Parameters &conditioned_on_parameters);
   */
  
  void subsample_specific_move(Ensemble* next_particles,
                               Ensemble* current_particles);
  
  /*
   void subsample_specific_move(Ensemble* next_particles,
   Ensemble* current_particles,
   const Parameters &conditioned_on_parameters);
   */
  
  //void pack(Ensemble &next_ensemble_kalman);
  //void unpack(Ensemble &next_ensemble_kalman);
  
  /**
   * @brief Copies the state of another SequentialEnsembleKalmanWorker into this object.
   *
   * @param another The SequentialEnsembleKalmanWorker instance to copy from.
   */
  void make_copy(const SequentialEnsembleKalmanWorker &another);
  
  /** @brief The log unnormalised incremental weights. */
  arma::colvec log_unnormalised_incremental_weights;
  
  /*
   void specific_move(EnsembleKalman* next_ensemble_kalman,
   const EnsembleKalman* current_ensemble_kalman);
   //void specific_simulate_and_weight(const Parameters &conditioned_on_parameters);
   
   void specific_move(EnsembleKalman* next_ensemble_kalman,
   const EnsembleKalman* current_ensemble_kalman,
   const Parameters &conditioned_on_parameters);
   
   void subsample_specific_simulate(EnsembleKalman* next_ensemble_kalman,
   const Parameters &conditioned_on_parameters);
   
   void subsample_specific_move(EnsembleKalman* next_ensemble_kalman,
   const EnsembleKalman* current_ensemble_kalman,
   const Parameters &conditioned_on_parameters);
   
   // Don't think we want to store this here, since found in ensemble_kalman (I think).
   //std::vector<Simulator*> simulate_priors;
   //std::vector<> simulate_for_likelihoods;
   
   // Not stored here.
   //EnsembleKalman* ensemble_kalman;
   
   arma::colvec log_unnormalised_incremental_weights;
   */
  
};

}

#endif
