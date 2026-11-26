#ifndef EXCHANGERATIOSTEP_H
#define EXCHANGERATIOSTEP_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <string>
#include <vector>
#include "ratio_step.h"
#include "unnormalised_normaliser_estimator.h"

namespace ilike
{

// One unnormalised distribution whose normalising constant the exchange move
// cancels: x ~ gamma(x | theta)/Z(theta) (dev/plan.md 3.22 U1)
struct ExchangeField
{
  std::string variable;
  UnnormalisedPotentialPtr potential = NULL;    // log gamma(field | parameters), up to a constant
  UnnormalisedExactFunctionPtr exact = NULL;    // an exact draw (NULL: sweeps)
  UnnormalisedSweepFunctionPtr sweep = NULL;    // one heat-bath sweep, in place
  size_t sweeps = 100;
  size_t replicates = 1;   // (replicates = R: one field per replicate in the target)
};

// The exchange algorithm's auxiliary draw (Murray, Ghahramani & MacKay 2006;
// dev/plan.md 3.22 U1): u ~ gamma(. | to)/Z(to) for each field, and the term
// gamma(u | from)/gamma(u | to), an unbiased estimate of Z(from)/Z(to), so
// that every normalising constant cancels in the acceptance ratio. The field
// may be data (observed) or a parameter (latent, read at its current value
// in `to`, which the sweeps start from). Exact when u is drawn exactly; with
// sweeps it is approximate (Everitt 2012).
class ExchangeRatioStep : public RatioStep
{
public:

  ExchangeRatioStep(const std::vector<ExchangeField> &fields_in,
                    Data* data_in,
                    const std::string &tempering_variable_in = "");

  RatioStep* duplicate() const;

  double draw(RandomNumberGenerator &rng,
              const Parameters &from,
              const Parameters &to) const;

protected:

  std::vector<ExchangeField> fields;
  // (dev/plan.md 3.22 U3b) under data tempering, the sequence's variable: the
  // target has R x fraction whole replicates, so as many fields are drawn
  std::string tempering_variable;

  // not stored here
  Data* data;

};

}

#endif
