#ifndef ADAPTATIONCRITERION_H
#define ADAPTATIONCRITERION_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <deque>
#include <memory>
#include <string>
#include <vector>

#include "recipe_value.h"

namespace ilike
{

class MCMCAdaptor;

// When an MCMC adaptor stops adapting (dev/plan.md 3.18, item 10; `until =`
// on an adaptor). Checked after each adaptation step; once it holds, the
// adaptor stops and the proposal keeps its parameters, so from then on the
// chain is a time-homogeneous Markov chain started where it stopped (valid
// whatever the criterion). Criteria:
// - iterations(n): after n adaptation steps;
// - stable(tolerance, window, watch, change): the watched quantities changed
//   by at most `tolerance` (relative or absolute) over the last `window`
//   steps (every value in the window against the current one);
// - acceptance(tolerance, window): the average acceptance probability over
//   the last `window` steps within `tolerance` of the adaptor's target;
// - all(…), any(…): of the criteria inside (each updated at every step).
// The state (the windows) belongs to one chain: copied with the adaptor.
class AdaptationCriterion
{
public:

  enum class Kind { iterations, stable, acceptance, all, any };

  AdaptationCriterion();
  AdaptationCriterion(const AdaptationCriterion &another);
  AdaptationCriterion &operator=(const AdaptationCriterion &another);
  AdaptationCriterion* duplicate() const;

  // from the recipe's node: {criterion, n | tolerance, window, watch, change | items}
  static AdaptationCriterion* from_recipe(const RecipeValue &node);

  // after adaptation step number `steps`, with that step's acceptance probability
  bool met(const MCMCAdaptor &adaptor, double acceptance_probability, size_t steps);
  // iterations(n): n; otherwise 0 (the length of the warm-up is not known in
  // advance; C220)
  size_t known_iterations() const { return (this->kind == Kind::iterations) ? this->n : 0; }

private:

  void make_copy(const AdaptationCriterion &another);

  Kind kind;
  size_t n;
  double tolerance;
  size_t window;
  std::vector<std::string> watch; // empty: the adaptor's own quantities
  bool relative;
  std::deque<std::vector<double>> history;
  std::deque<double> acceptances;
  std::vector<std::unique_ptr<AdaptationCriterion>> items;
};

}

#endif
