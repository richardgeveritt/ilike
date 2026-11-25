#ifndef CALLER_TEMPERING_H
#define CALLER_TEMPERING_H

#include <cmath>
#include <cstddef>

// (dev/plan.md 3.21 (c2'), (c2''), (c2''')) A filter inside a factor that the
// caller's SMC takes in one observation at a time (data tempering,
// data(fraction, ...)): the rules shared by the particle, Kalman and ensemble
// Kalman filters, so that they cannot drift apart.

namespace ilike
{

// The caller's target: n fraction = k + r, k whole observations and a share
// r in [0, 1) of the next. Without anneal, k = ceil(n fraction) and r = 0, as
// data(..., anneal = false) adds whole data points to factors of model
// statements; with anneal, k = floor(n fraction), as plates do. A fraction
// written as j/n gives k = j, r = 0 in both (rounding allowed for).
inline void caller_tempering_point(double fraction,
                                   std::size_t number_of_observations,
                                   bool anneal,
                                   std::size_t &k,
                                   double &r)
{
  double n = double(number_of_observations);
  double x = n*fraction;
  double tiny = 1e-9*n;
  double whole;
  r = 0.0;
  if (!anneal)
    whole = std::ceil(x - tiny);
  else
  {
    whole = std::floor(x + tiny);
    if (x - whole > tiny)
      r = x - whole;
  }
  if (whole < 0.0) { whole = 0.0; r = 0.0; }
  if (whole >= n) { whole = n; r = 0.0; }
  k = std::size_t(whole);
}

// What a trial at the caller's target (k, r) is, for an output that has
// committed `committed` observations
enum class CallerTrial
{
  committed,   // the committed observations only (k == committed, r == 0; or k < committed)
  next,        // the next observation, whole (k == committed + 1, r == 0) or a share r (k == committed, r > 0)
  too_far      // more than one observation beyond the committed ones: refused
};

inline CallerTrial caller_trial(std::size_t committed, std::size_t k, double r)
{
  if ((k == committed) && (r == 0.0))
    return CallerTrial::committed;
  if (((k == committed) && (r > 0.0)) || ((k == committed + 1) && (r == 0.0)))
    return CallerTrial::next;
  if (k < committed)
    return CallerTrial::committed; // (the caller's fraction never decreases; kept as before)
  return CallerTrial::too_far;
}

// (dev/plan.md 3.21 (c2''), F2; C267) the exchange step's new size for a
// filter with `current` particles (or ensemble members), or 0 for no exchange:
// only when the caller's acceptance rate is below the threshold, the size is
// below the limit (max 0: none), and the caller's target is at a whole
// observation (r == 0: between two, the targets for different sizes have
// different marginals, C258)
inline std::size_t caller_exchange_size(std::size_t current,
                                        double acceptance_rate,
                                        double below,
                                        double multiply,
                                        std::size_t maximum,
                                        double r)
{
  if ((below < 0.0) || (acceptance_rate >= below) || (r > 0.0))
    return 0;
  if ((maximum > 0) && (current >= maximum))
    return 0;
  std::size_t new_size = std::size_t(std::ceil(multiply*double(current)));
  if (new_size <= current)
    new_size = current + 1;
  if ((maximum > 0) && (new_size > maximum))
    new_size = maximum;
  return new_size;
}

// the message for CallerTrial::too_far
inline const char* caller_trial_too_far_message()
{
  return "SMC^2: data tempering of a filter takes in at most one observation per step: each whole observation "
         "(a fraction j/n) must be a point of the sequence's schedule.";
}

}

#endif
