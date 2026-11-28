#ifndef HILBERTSORT_H
#define HILBERTSORT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <cstdint>

namespace ilike
{

// (dev/plan.md 3.22 U9 R3c, C307) the order in which to resample the
// particles of a correlated particle filter (Deligiannidis, Doucet & Pitt
// 2018): points close in space close in the order, so that a small change in
// the particles (from a small change in the parameters or the random numbers)
// changes the ancestors little. points is n x k (one row per particle). k = 1:
// by value. k > 1: each coordinate standardised (mean and sd over the
// points) and mapped to (0, 1) by the logistic function (Gerber & Chopin's
// SQMC), then ordered by the index of its cell along the Hilbert curve, with
// min(12, 64/k) bits per coordinate (Skilling 2004, "Programming the Hilbert curve").
// Ties keep the particles' order.
std::vector<size_t> hilbert_order(const arma::mat &points);

// the Hilbert index of a point of [0, 2^bits)^n given as integers (n * bits
// <= 64); x is overwritten
uint64_t hilbert_index(uint64_t* x, size_t n, unsigned bits);

}

#endif
