#ifndef PRINCIPALCOMPONENTSINFO_H
#define PRINCIPALCOMPONENTSINFO_H

#include <RcppArmadillo.h>

namespace ilike
{

// The state of a principal-components random walk (Andrieu & Thoms 2008,
// Algorithm 8): m orthonormal directions w_l (the columns of directions),
// the variances rho_l of the target along them, and a scale ell_l for each.
// A move picks a direction l and adds N(0, ell_l rho_l) w_l.
struct PrincipalComponentsInfo
{
  arma::mat directions;
  arma::colvec eigenvalues;
  arma::colvec scales;
};

}

#endif
