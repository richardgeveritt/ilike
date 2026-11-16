#ifndef SAMPLE_STORE_H
#define SAMPLE_STORE_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <string>
#include <boost/any.hpp>
#include "parameters.h"

namespace ilike
{

// Compact, columnar storage of a sequence of points (dev/plan.md 2.2, 2.5
// "Stage 2 - concrete design"): what is kept of the points an algorithm no
// longer works on (e.g. the iterations of an MCMC chain), instead of one
// Particle each.
//
// The layout is fixed by the first point appended: every matrix (vector)
// variable gets a block of columns, filled with the variable's elements in
// the order of Parameters::get_rowvec (Armadillo's as_row), and its shape is
// kept so that parameters_at can rebuild it. Points are stored as the columns
// of one contiguous matrix that grows by doubling. Non-matrix ("any")
// variables are kept as boost::any objects, so that nothing assumes every
// variable is numeric. Algorithms may also declare extra numeric columns per
// point (e.g. a log weight, or one value per factor).
class SampleStore
{
public:

  SampleStore();
  virtual ~SampleStore();

  SampleStore(const SampleStore &another);
  SampleStore& operator=(const SampleStore &another);
  SampleStore(SampleStore &&another);
  SampleStore& operator=(SampleStore &&another);

  // Declare number_of_extra_columns per-point extra values. Only before the
  // first point is appended.
  void set_number_of_extra_columns(size_t number_of_extra_columns_in);
  size_t get_number_of_extra_columns() const;

  // True if the point has exactly the variables, shapes and kinds of the
  // stored layout (always true while the store is empty).
  bool compatible(const Parameters &parameters) const;

  // Append a point (and its extra values, if extra columns were declared).
  // Stops if the point is not compatible with the layout.
  void append(const Parameters &parameters,
              const double *extra_values = NULL);

  size_t size() const;
  bool empty() const;
  void reserve(size_t number_of_points);
  void clear();

  // The points as a (points x columns) matrix, with the columns of the given
  // variables in the given order: the same as stacking
  // Parameters::get_rowvec(variables) for every point.
  arma::mat matrix_of(const std::vector<std::string> &variables) const;

  // Copy the given variables of every point into columns first_column, ...
  // of destination (one column per point): the transpose of matrix_of.
  void fill_transposed(const std::vector<std::string> &variables,
                       arma::mat &destination,
                       size_t first_column) const;

  // Reserve room for number_of_points points once the layout is known (the
  // request is remembered until the first point is appended).
  void reserve_after_first(size_t number_of_points);

  // The stored point at index i, rebuilt as a Parameters object.
  Parameters parameters_at(size_t i) const;

  // Extra value j of point i.
  double extra_at(size_t i, size_t j) const;

  // The extra values of point i (contiguous).
  const double* extras_of(size_t i) const;

  // Bytes held by the numeric storage (for benchmarks).
  size_t memory_bytes() const;

protected:

  void setup_layout(const Parameters &parameters);
  void grow(size_t minimum_capacity);

  bool layout_set;

  // matrix variables: name, first column, number of elements, shape
  std::vector<std::string> vector_names;
  std::vector<size_t> vector_offsets;
  std::vector<arma::uword> vector_n_rows;
  std::vector<arma::uword> vector_n_cols;
  size_t total_dimension;

  // any variables: name and values
  std::vector<std::string> any_names;
  std::vector<std::vector<boost::any>> any_values;

  size_t number_of_extra_columns;

  size_t number_of_points;
  size_t pending_reserve = 0;

  // Where each stored variable was found in the last Parameters appended
  // (its position among the variables), so that append looks there first
  // instead of searching by name (dev/plan.md, Stage 3 item 7, in the store
  // only). Checked by name on every use; a miss searches and updates it.
  std::vector<size_t> vector_positions;
  const arma::mat& find_vector_variable(const Parameters &parameters, size_t k);

  // (total_dimension x capacity) and (number_of_extra_columns x capacity):
  // one column per point
  arma::mat values;
  arma::mat extras;

};

}

#endif
