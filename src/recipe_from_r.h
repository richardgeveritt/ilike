#ifndef RECIPE_FROM_R_H
#define RECIPE_FROM_R_H

#include <RcppArmadillo.h>
#include "recipe_value.h"

namespace ilike
{

// R binding: an R value (a recipe made by compile(), a list of run-time
// constants, a method list) as a RecipeValue (dev/plan.md 3.14, I1).
RecipeValue recipe_from_R(SEXP x);

// and back (functions as "language::name", values of other kinds as NULL)
SEXP recipe_to_R(const RecipeValue &x);

}

#endif
