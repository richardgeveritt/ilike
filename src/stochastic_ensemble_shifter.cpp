#include "stochastic_ensemble_shifter.h"
#include "distributions.h"
#include "ensemble_factor_variables.h"

namespace ilike {
StochasticEnsembleShifter::StochasticEnsembleShifter() : EnsembleShifter() {}

StochasticEnsembleShifter::~StochasticEnsembleShifter() {}

// Copy constructor for the StochasticEnsembleShifter class.
StochasticEnsembleShifter::StochasticEnsembleShifter(
    const StochasticEnsembleShifter &another)
    : EnsembleShifter(another) {
  this->make_copy(another);
}

void StochasticEnsembleShifter::operator=(
    const StochasticEnsembleShifter &another) {
  if (this == &another) { // if a==a
    return;
  }

  EnsembleShifter::operator=(another);
  this->make_copy(another);
}

EnsembleShifter *StochasticEnsembleShifter::duplicate() const {
  return (new StochasticEnsembleShifter(*this));
}

void StochasticEnsembleShifter::make_copy(
    const StochasticEnsembleShifter &another) {}

void StochasticEnsembleShifter::setup(Ensemble *ensemble,
                                      double inverse_incremental_temperature) {}

void StochasticEnsembleShifter::shift(
    const std::vector<arma::colvec> &shift_terms,
    arma::colvec &position, const std::vector<arma::colvec *> &measurements,
    const std::vector<arma::mat> &kalman_gains,
    const std::vector<arma::colvec> &myys,
    double inverse_incremental_temperature) const {
  for (size_t j = 0; j < kalman_gains.size(); ++j) {
    position = position + kalman_gains[j] * (*measurements[j] - shift_terms[j]);
  }
}
} // namespace ilike
