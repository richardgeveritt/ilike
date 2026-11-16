#ifndef SHARED_POINT_WEIGHTS_H
#define SHARED_POINT_WEIGHTS_H

#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include "parameters.h"

namespace ilike
{

// Points of several particles, each of which may hold several points (e.g.
// the states of a chain), with particle i's normalised weight shared equally
// between its k_i points: each gets log weight normalised_log_weights[i] -
// log(k_i). The weights of all the points then still sum to one, as mean_wt
// and cov_wt require. Particles without points are skipped (their weight is
// lost, so the weights then sum to less than one; no MoveOutput has none).
inline void shared_point_weights(const std::vector< std::vector<Parameters> > &points_per_particle,
                                 const arma::colvec &normalised_log_weights,
                                 std::vector<Parameters> &all_points,
                                 arma::colvec &all_normalised_log_weights)
{
  if (normalised_log_weights.n_elem!=points_per_particle.size())
    Rcpp::stop("shared_point_weights: one weight per particle needed.");
  size_t total = 0;
  for (size_t i=0; i<points_per_particle.size(); ++i)
    total = total + points_per_particle[i].size();
  all_points.clear();
  all_points.reserve(total);
  all_normalised_log_weights.set_size(total);
  size_t counter = 0;
  for (size_t i=0; i<points_per_particle.size(); ++i)
  {
    size_t k = points_per_particle[i].size();
    if (k==0)
      continue;
    double share = normalised_log_weights[i] - std::log(double(k));
    for (size_t j=0; j<k; ++j)
    {
      all_points.push_back(points_per_particle[i][j]);
      all_normalised_log_weights[counter] = share;
      counter = counter + 1;
    }
  }
}

}

#endif
