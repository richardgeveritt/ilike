#ifndef ITERATIVEENSEMBLEKALMANINVERSION_H
#define ITERATIVEENSEMBLEKALMANINVERSION_H

#include "ensemble_kalman.h"
#include "ensemble_sequencer.h"

class EnsembleKalmanOutput;
class MCMC;
class MoveOutput;

class IterativeEnsembleKalmanInversion : public EnsembleKalman
{
public:

  IterativeEnsembleKalmanInversion();
  
  IterativeEnsembleKalmanInversion(const IterativeEnsembleKalmanInversion &another);
  
  virtual ~IterativeEnsembleKalmanInversion(void);

  void operator=(const IterativeEnsembleKalmanInversion &another);
  EnsembleKalman* ensemble_kalman_duplicate() const;
  LikelihoodEstimator* duplicate() const;
  
  MoveOutput* move(RandomNumberGenerator &rng,
                   Particle &particle);
  
  MoveOutput* move(RandomNumberGenerator &rng,
                   Particle &particle,
                   const Parameters &conditioned_on_parameters);
  
  MoveOutput* subsample_move(RandomNumberGenerator &rng,
                             Particle &particle,
                             const Parameters &conditioned_on_parameters);
  
  //void weight_for_adapting_sequence(Particles &current_particles);
  
  //void weight_for_adapting_sequence(Particles &current_particles,
  //                                  const Parameters &conditioned_on_parameters);
  
  //void subsample_weight_for_adapting_sequence(Particles &current_particles,
  //                                            const Parameters &conditioned_on_parameters);
  
  //void weight_for_adapting_sequence(Ensemble &current_particles,
  //                                  double incremental_temperature);

protected:
  
  EnsembleKalmanOutput* specific_run();
  
  EnsembleKalmanOutput* ensemble_kalman_initialise();
  void ensemble_kalman_simulate(EnsembleKalmanOutput* simulation);
  void ensemble_kalman_evaluate(EnsembleKalmanOutput* simulation);
  void ensemble_kalman_evaluate_smcfixed_part(EnsembleKalmanOutput* simulation);
  void ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput* simulation);
  
  EnsembleKalmanOutput* specific_run(const Parameters &parameters);
  
  EnsembleKalmanOutput* ensemble_kalman_initialise(const Parameters &parameters);
  void ensemble_kalman_simulate(EnsembleKalmanOutput* simulation,
                                const Parameters &conditioned_on_parameters);
  
  void ensemble_kalman_evaluate(EnsembleKalmanOutput* simulation,
                                const Parameters &conditioned_on_parameters);
  void ensemble_kalman_evaluate_smcfixed_part(EnsembleKalmanOutput* simulation,
                                              const Parameters &conditioned_on_parameters);
  void ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput* simulation,
                                                                const Parameters &conditioned_on_parameters);
  
  void ensemble_kalman_subsample_simulate(EnsembleKalmanOutput* simulation,
                                          const Parameters &conditioned_on_parameters);
  
  void ensemble_kalman_subsample_evaluate(EnsembleKalmanOutput* simulation,
                                          const Parameters &conditioned_on_parameters);
  void ensemble_kalman_subsample_evaluate_smcfixed_part(EnsembleKalmanOutput* simulation,
                                                        const Parameters &conditioned_on_parameters);
  void ensemble_kalman_subsample_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput* simulation,
                                                                         const Parameters &conditioned_on_parameters);

  //void smc_update(EnsembleKalmanOutput* current_state);

  void make_copy(const IterativeEnsembleKalmanInversion &another);
  
  // Stored here.
  MCMC* mcmc;
  
  EnsembleSequencer sequencer;
  bool sequencer_limit_is_fixed;
  
  // stored here
  Index* index;

  //void smc_step();

  //void weight_update();
};

#endif
