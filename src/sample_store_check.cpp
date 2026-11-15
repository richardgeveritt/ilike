#include <RcppArmadillo.h>
#include "sample_store.h"
using namespace Rcpp;
using namespace ilike;

// Internal check of SampleStore (tests/testthat/test-sample-store.R): builds
// a store from points with a matrix variable, a scalar and a non-matrix
// variable, and returns what the store gives next to what the existing
// Parameters functions give, for the test to compare.
// [[Rcpp::export]]
List sample_store_check()
{
  std::vector<Parameters> points;
  for (int i = 0; i < 5; ++i)
  {
    Parameters p;
    arma::mat theta(2, 3);
    for (arma::uword k = 0; k < theta.n_elem; ++k) theta[k] = 10.0*i + k;
    p["theta"] = theta;
    p["sigma"] = arma::mat(1, 1, arma::fill::value(-1.0*i));
    p("label") = std::string("point") + std::to_string(i);
    points.push_back(p);
  }

  SampleStore store;
  store.set_number_of_extra_columns(2);
  store.reserve_after_first(2);   // smaller than needed: must still grow
  for (int i = 0; i < 5; ++i)
  {
    double extra[2] = {double(i), 100.0 + i};
    store.append(points[i], extra);
  }

  std::vector<std::string> variables = {"sigma", "theta"};
  arma::mat expected(5, 7);
  for (int i = 0; i < 5; ++i)
    expected.row(i) = points[i].get_rowvec(variables);

  arma::mat transposed(7, 8, arma::fill::zeros);
  store.fill_transposed(variables, transposed, 2);

  Parameters back = store.parameters_at(3);
  arma::mat theta3 = back["theta"];
  arma::mat sigma3 = back["sigma"];
  std::string label3 = boost::any_cast<std::string>(back("label"));

  Parameters wrong_size;
  wrong_size["theta"] = arma::mat(3, 2, arma::fill::zeros);
  wrong_size["sigma"] = arma::mat(1, 1, arma::fill::zeros);
  wrong_size("label") = std::string("x");
  Parameters missing;
  missing["theta"] = arma::mat(2, 3, arma::fill::zeros);
  missing("label") = std::string("x");

  arma::mat extras(5, 2);
  for (int i = 0; i < 5; ++i) for (int j = 0; j < 2; ++j) extras(i, j) = store.extra_at(i, j);

  SampleStore copy = store;
  copy.clear();

  return List::create(_["size"] = store.size(),
                      _["matrix"] = store.matrix_of(variables),
                      _["expected"] = expected,
                      _["transposed"] = transposed,
                      _["theta3"] = theta3,
                      _["expected_theta3"] = arma::mat(points[3]["theta"]),
                      _["sigma3"] = sigma3,
                      _["label3"] = label3,
                      _["extras"] = extras,
                      _["compatible_same"] = store.compatible(points[0]),
                      _["compatible_wrong_size"] = store.compatible(wrong_size),
                      _["compatible_missing"] = store.compatible(missing),
                      _["copy_cleared_size"] = copy.size(),
                      _["original_size_after_copy_cleared"] = store.size());
}
