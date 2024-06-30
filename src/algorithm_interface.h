//
//  algorithm_interface.h
//  ilike_cpp
//
//  Created by Richard Everitt on 06/01/2024.
//

#ifndef algorithm_interface_h
#define algorithm_interface_h

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "distributions.h"
#include "parameters.h"

class KalmanFilter;
class EnsembleKalmanFilter;
class ParticleFilter;
class EnsembleShifter;

KalmanFilter* get_kalman_filter(Data* the_data,
                                const List &model,
                                const List &parameters,
                                const List &algorithm_parameter_list,
                                size_t kf_iterations_to_store,
                                bool write_to_file_at_each_iteration,
                                const std::string &results_name_in);

EnsembleKalmanFilter* get_ensemble_kalman_filter(RandomNumberGenerator* rng,
                                                 Data* the_data,
                                                 const List &model,
                                                 const List &parameters,
                                                 const List &algorithm_parameter_list,
                                                 size_t number_of_ensemble_members,
                                                 EnsembleShifter* shifter,
                                                 size_t enk_iterations_to_store,
                                                 bool write_to_file_at_each_iteration,
                                                 bool parallel_in,
                                                 size_t grain_size_in,
                                                 const String &results_name_in,
                                                 size_t* seed,
                                                 std::vector<Data> &data_created_in_get_measurement_covariance_estimators);

ParticleFilter* get_particle_filter(RandomNumberGenerator* rng,
                                    Data* the_data,
                                    const List &model,
                                    const List &parameters,
                                    const List &algorithm_parameter_list,
                                    size_t number_of_particles,
                                    size_t smc_iterations_to_store,
                                    bool write_to_file_at_each_iteration,
                                    bool parallel_in,
                                    size_t grain_size_in,
                                    const String &results_name_in,
                                    size_t* seed,
                                    std::vector<Data> &data_created_in_get_likelihood_estimators,
                                    std::vector<Data> &data_created_in_get_measurement_covariance_estimators);

#endif /* algorithm_interface_h */
