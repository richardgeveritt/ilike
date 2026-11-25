#include "exact_kalman_updater.h"
#include "utils.h"
#include "kalman_filter_output.h"

namespace ilike
{

// (C276) exactly the same matrix (sizes and every value)
static bool same_matrix(const arma::mat &a, const arma::mat &b)
{
  return (a.n_rows == b.n_rows) && (a.n_cols == b.n_cols) &&
    std::equal(a.memptr(), a.memptr() + a.n_elem, b.memptr());
}

// (C274) A = (A + A')/2 without a temporary: the same numbers as
// 0.5*(A + A.t())
static void symmetrise_in_place(arma::mat &A)
{
  for (arma::uword j = 0; j < A.n_cols; ++j)
    for (arma::uword i = j + 1; i < A.n_rows; ++i)
    {
      double v = 0.5*(A(i, j) + A(j, i));
      A(i, j) = v;
      A(j, i) = v;
    }
}
ExactKalmanUpdater::ExactKalmanUpdater()
:KalmanUpdater()
{
}

ExactKalmanUpdater::ExactKalmanUpdater(const std::string &state_variable_in,
                                       const std::string &measurement_variable_in,
                                       const arma::mat &measurement_matrix_in,
                                       const arma::mat &measurement_noise_in)
:KalmanUpdater(state_variable_in,measurement_variable_in)
{
  this->measurement_matrix = measurement_matrix_in;
  this->measurement_noise = measurement_noise_in;
  this->set_using_parameters = false;
}

ExactKalmanUpdater::ExactKalmanUpdater(const std::string &state_variable_in,
                                       const std::string &measurement_variable_in,
                                       GetMatrixPtr measurement_matrix_function_in,
                                       GetMatrixPtr measurement_noise_function_in)
:KalmanUpdater(state_variable_in,measurement_variable_in)
{
  this->measurement_matrix_function = measurement_matrix_function_in;
  this->measurement_noise_function = measurement_noise_function_in;
  this->set_using_parameters = true;
}

ExactKalmanUpdater::~ExactKalmanUpdater()
{
  
}

ExactKalmanUpdater::ExactKalmanUpdater(const ExactKalmanUpdater &another)
:KalmanUpdater(another)
{
  this->make_copy(another);
}

void ExactKalmanUpdater::operator=(const ExactKalmanUpdater &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  KalmanUpdater::operator=(another);
  this->make_copy(another);
}

KalmanUpdater* ExactKalmanUpdater::duplicate() const
{
  return( new ExactKalmanUpdater(*this));
}

void ExactKalmanUpdater::make_copy(const ExactKalmanUpdater &another)
{
  this->measurement_matrix = another.measurement_matrix;
  this->measurement_noise = another.measurement_noise;
  this->measurement_matrix_function = another.measurement_matrix_function;
  this->measurement_noise_function = another.measurement_noise_function;
  this->tracks_matrix_changes = another.tracks_matrix_changes; // (C276)
  this->changed_at_last_set = another.changed_at_last_set;
  this->frozen_gain = another.frozen_gain;
  this->frozen_L = another.frozen_L;
  this->frozen_log_constant = another.frozen_log_constant;
}

void ExactKalmanUpdater::update(KalmanFilterOutput* current_state,
                                const arma::colvec &current_measurement)
{
  // (C271) from the current posterior, which the filter sets to the
  // prediction before the first update: a second measurement factor updates
  // the first one's posterior (each used to start again from the stored
  // prediction). The matrices are those set_parameters evaluated (they were
  // evaluated again here).
  const arma::colvec &predicted_mean = current_state->get_current_posterior_mean();
  const arma::mat &predicted_covariance = current_state->get_current_posterior_covariance();
  // (C272) the update written around G = H P, with one Cholesky factor of
  // the innovation covariance S = G H' + R = L L':
  //   W = L^{-1} G,  z = L^{-1} v  (v = y - H mean, the innovation)
  //   mean+ = mean + G' S^{-1} v = mean + W' z
  //   P+    = P - G' S^{-1} G   = P - W' W
  //   log N(y; H mean, S) = -(d/2) log(2 pi) - sum(log diag L) - z'z/2
  // This is the same Kalman update as K = P H' S^{-1}, P+ = (I - K H) P,
  // without the m x m x m product (I - K H) P and with one factorisation of
  // S instead of three (a general inverse, then log_det_sympd and inv_sympd
  // in dmvnorm). S and P+ are made exactly symmetric (the products leave
  // them symmetric only up to rounding).
  // (C274) the temporaries are members, so their memory is reused between
  // steps; with one measurement component (d = 1) S is a number, and L its
  // square root, so no factorisation or triangular solve is needed (the same
  // arithmetic: chol of a 1 x 1 matrix is the square root, its triangular
  // solves divisions).
  const arma::mat &H = this->measurement_matrix;
  double incremental_log_llhd;
  if ((H.n_rows > 1) && this->measurement_noise.is_diagmat())
  {
    incremental_log_llhd = this->update_one_component_at_a_time(predicted_mean, predicted_covariance, current_measurement);
    current_state->set_current_posterior_statistics(this->new_mean, this->new_covariance);
    current_state->add_log_normalising_constant_ratio(incremental_log_llhd);
    current_state->log_likelihood = current_state->log_likelihood + incremental_log_llhd;
    return;
  }
  this->innovation = current_measurement - H*predicted_mean;
  this->G = H*predicted_covariance;
  if (H.n_rows == 1)
  {
    double s = arma::dot(this->G, H) + this->measurement_noise(0, 0);
    if (!(s > 0.0))
      Rcpp::stop("kf: the variance of the predicted measurement of '" + this->measurement_variable +
                 "' is not positive (check the data model's covariance).");
    double l = std::sqrt(s);
    this->W = this->G/l;
    double z0 = this->innovation(0)/l;
    this->new_mean = predicted_mean + z0*this->W.t();
    incremental_log_llhd = -0.5*log(2.0*M_PI) - log(l) - 0.5*z0*z0;
  }
  else
  {
    this->S = this->G*H.t() + this->measurement_noise;
    symmetrise_in_place(this->S);
    if (!arma::chol(this->L, this->S, "lower"))
      Rcpp::stop("kf: the covariance of the predicted measurement of '" + this->measurement_variable +
                 "' is not positive definite (check the data model's covariance).");
    arma::solve(this->W, arma::trimatl(this->L), this->G);
    arma::solve(this->z, arma::trimatl(this->L), this->innovation);
    this->new_mean = predicted_mean + this->W.t()*this->z;
    incremental_log_llhd = -0.5*double(this->innovation.n_elem)*log(2.0*M_PI) -
      arma::accu(arma::log(this->L.diag())) - 0.5*arma::dot(this->z, this->z);
  }
  this->new_covariance = predicted_covariance - this->W.t()*this->W;
  symmetrise_in_place(this->new_covariance);

  current_state->set_current_posterior_statistics(this->new_mean, this->new_covariance);
  current_state->add_log_normalising_constant_ratio(incremental_log_llhd);
  current_state->log_likelihood = current_state->log_likelihood + incremental_log_llhd;
}

// (C275) With a diagonal measurement covariance R the components of y are
// independent given the state, so conditioning on y is conditioning on y_1,
// then y_2 given y_1, and so on (the univariate treatment of Koopman & Durbin,
// 2000): for each component j, with h the j-th row of H,
//   g = h P,  s = g h' + R_jj,  v = y_j - h mean,
//   mean += g' v / s,  P -= g' g / s,  log-density += log N(v; 0, s).
// The same posterior and log-likelihood as the update with all components
// together (up to rounding), without forming or factorising the d x d
// matrix S: about 2 m^2 d operations, against about 2 m^2 d + 2 m d^2 + d^3/3.
// P -= g' g / s is done in place, and P made exactly symmetric at the end.
double ExactKalmanUpdater::update_one_component_at_a_time(const arma::colvec &predicted_mean,
                                                          const arma::mat &predicted_covariance,
                                                          const arma::colvec &current_measurement)
{
  const arma::mat &H = this->measurement_matrix;
  this->new_mean = predicted_mean;
  this->new_covariance = predicted_covariance;
  symmetrise_in_place(this->new_covariance);
  arma::uword m = this->new_mean.n_elem;
  double log_llhd = 0.0;
  for (arma::uword j = 0; j < H.n_rows; ++j)
  {
    this->g = H.row(j)*this->new_covariance;
    double s = arma::dot(this->g, H.row(j)) + this->measurement_noise(j, j);
    if (!(s > 0.0))
      Rcpp::stop("kf: the variance of component " + std::to_string(j + 1) + " of the predicted measurement of '" +
                 this->measurement_variable + "' is not positive (check the data model's covariance).");
    double v = current_measurement(j) - arma::dot(H.row(j), this->new_mean);
    double v_over_s = v/s;
    for (arma::uword r = 0; r < m; ++r)
      this->new_mean(r) += this->g(r)*v_over_s;
    for (arma::uword c = 0; c < m; ++c)
    {
      double gc_over_s = this->g(c)/s;
      for (arma::uword r = 0; r < m; ++r)
        this->new_covariance(r, c) -= this->g(r)*gc_over_s;
    }
    log_llhd += -0.5*(log(2.0*M_PI) + log(s) + v*v_over_s);
  }
  symmetrise_in_place(this->new_covariance);
  return log_llhd;
}

void ExactKalmanUpdater::set_parameters(const Parameters &conditioned_on_parameters_in)
{
  if (this->set_using_parameters)
  {
    if (this->tracks_matrix_changes)
    {
      arma::mat H = this->measurement_matrix_function(conditioned_on_parameters_in);
      arma::mat R = this->measurement_noise_function(conditioned_on_parameters_in);
      this->changed_at_last_set = !same_matrix(H, this->measurement_matrix) || !same_matrix(R, this->measurement_noise);
      this->measurement_matrix = std::move(H);
      this->measurement_noise = std::move(R);
    }
    else
    {
      this->measurement_matrix = this->measurement_matrix_function(conditioned_on_parameters_in);
      this->measurement_noise = this->measurement_noise_function(conditioned_on_parameters_in);
    }
  }
}

// (C276) for the covariance P entering this update (once the filter is in
// steady state, the same at every step): S = H P H' + R = L L', and the gain
// K = P H' S^{-1} (as K' = L'^{-1} L^{-1} H P), from the full formula whichever
// path update takes
void ExactKalmanUpdater::freeze_gain(const arma::mat &entering_covariance)
{
  const arma::mat &H = this->measurement_matrix;
  arma::mat G = H*entering_covariance;
  arma::mat S = G*H.t() + this->measurement_noise;
  symmetrise_in_place(S);
  if (!arma::chol(this->frozen_L, S, "lower"))
    Rcpp::stop("kf: the covariance of the predicted measurement of '" + this->measurement_variable +
               "' is not positive definite (check the data model's covariance).");
  arma::mat W = arma::solve(arma::trimatl(this->frozen_L), G);
  this->frozen_gain = arma::solve(arma::trimatu(this->frozen_L.t()), W).t();
  this->frozen_log_constant = -0.5*double(H.n_rows)*log(2.0*M_PI) - arma::accu(arma::log(this->frozen_L.diag()));
}

// (C276) in steady state: mean += K v, log-density log N(v; 0, S) with the
// frozen K and L (v = y - H mean); the filter keeps the covariance
void ExactKalmanUpdater::update_mean(KalmanFilterOutput* current_state,
                                     const arma::colvec &current_measurement)
{
  const arma::colvec &mean = current_state->get_current_posterior_mean();
  this->innovation = current_measurement - this->measurement_matrix*mean;
  arma::solve(this->z, arma::trimatl(this->frozen_L), this->innovation);
  this->new_mean = mean + this->frozen_gain*this->innovation;
  double incremental_log_llhd = this->frozen_log_constant - 0.5*arma::dot(this->z, this->z);
  current_state->set_current_posterior_mean(this->new_mean);
  current_state->add_log_normalising_constant_ratio(incremental_log_llhd);
  current_state->log_likelihood = current_state->log_likelihood + incremental_log_llhd;
}
}
