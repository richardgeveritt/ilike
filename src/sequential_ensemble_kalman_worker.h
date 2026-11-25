#ifndef SEQUENTIALENSEMBLEKALMANWORKER_H
#define SEQUENTIALENSEMBLEKALMANWORKER_H

#include "ensemble_kalman_worker.h"

#include <vector>

#include "ensemble_kalman.h"

namespace ilike
{

class SequentialEnsembleKalmanWorker : public EnsembleKalmanWorker
{
public:
  
  SequentialEnsembleKalmanWorker();
  virtual ~SequentialEnsembleKalmanWorker();
  
  SequentialEnsembleKalmanWorker(EnsembleKalman* the_enk_in);
  
  SequentialEnsembleKalmanWorker(const SequentialEnsembleKalmanWorker &another);
  void operator=(const SequentialEnsembleKalmanWorker &another);
  EnsembleKalmanWorker* duplicate() const;
  
  void shift(Ensemble* ensemble,
             double inverse_incremental_temperature);
  
  void pack(Ensemble* ensemble);
  void unpack(Ensemble* ensemble);
  void unpack_with_predicted(Ensemble* ensemble);
  
  void weight(Ensemble* ensemble,
              const Index* index,
              double incremental_temperature);
  void subsample_weight(Ensemble* ensemble,
                        const Index* index,
                        double incremental_temperature);
  
  arma::colvec get_unnormalised_log_incremental_weights() const;

protected:
  
  void specific_simulate(Ensemble* next_ensemble_kalman,
                         const Index* index);
  
  void specific_simulate(Ensemble* next_ensemble_kalman,
                         const Index* index,
                         const Parameters &conditioned_on_parameters);
  
  void specific_move(Ensemble* next_particles,
                     Ensemble* current_particles);

  void subsample_specific_move(Ensemble* next_particles,
                               Ensemble* current_particles);

  void make_copy(const SequentialEnsembleKalmanWorker &another);
  
  arma::colvec log_unnormalised_incremental_weights;
public:
  // (C267) after the filter's number of ensemble members changed
  void resize_for_members() { this->log_unnormalised_incremental_weights.set_size(this->get_number_of_ensemble_members()); }

};

}

#endif
