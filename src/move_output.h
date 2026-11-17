#ifndef MOVEOUTPUT_H
#define MOVEOUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <iostream>
#include "particle.h"
#include "ilike_header.h"

namespace ilike
{

class Ensemble;
class Particles;

class MoveOutput
{
  
public:
  
  MoveOutput();
  virtual ~MoveOutput();
  
  MoveOutput(const MoveOutput &another);
  
  void operator=(const MoveOutput &another);
  virtual MoveOutput* duplicate() const=0;
  
  virtual Particle& back()=0;
  virtual Particle back() const=0;
  
  virtual std::vector<Parameters> get_vector_of_parameters() const=0;
  
  virtual void write_vector_points(const std::vector<std::string> &variables,
                                   std::ofstream &file_stream,
                                   std::shared_ptr<Transform> transform) const=0;
  virtual void write_any_points(const std::vector<std::string> &variables,
                                std::ofstream &file_stream) const=0;

  // Returns all vector points as a matrix (one row per MCMC step / particle).
  virtual arma::mat get_matrix_of_vector_points(const std::vector<std::string> &variables,
                                                std::shared_ptr<Transform> transform) const=0;

  // Write this output's points into columns first_column, first_column + 1,
  // ... of destination (one column per point, the variables' values down
  // the column): the transpose of get_matrix_of_vector_points(variables,
  // NULL). Used to write all points of a generation without intermediate
  // copies. Returns the number of columns written.
  virtual size_t fill_vector_points_transposed(const std::vector<std::string> &variables,
                                               arma::mat &destination,
                                               size_t first_column) const
  {
    arma::mat m = this->get_matrix_of_vector_points(variables, NULL);
    if (m.n_rows > 0)
      destination.cols(first_column, first_column + m.n_rows - 1) = m.t();
    return m.n_rows;
  }
  
  // The adaptation trace (one row per stored point, columns as named by
  // adaptation_names) and the state at the end of the chain (one row), for
  // outputs of moves that adapt (dev/plan.md 3.2); none by default.
  virtual size_t adaptation_trace_columns() const { return 0; }
  virtual size_t adaptation_final_columns() const { return 0; }
  virtual size_t fill_adaptation_trace_transposed(arma::mat &destination,
                                                  size_t first_column) const { return 0; }
  virtual void fill_adaptation_final(arma::mat &destination,
                                     size_t column) const {}
  virtual void adaptation_names(std::vector<std::string> &trace_names,
                                std::vector<std::string> &final_names) const {}

  virtual void write_factors(const std::string &directory_name,
                             const std::string &index) const=0;
  virtual void write_ensemble_factors(const std::string &directory_name,
                                      const std::string &index) const=0;
  
  virtual size_t length() const=0;

  // Make the factor variables of every particle held here belong to
  // `factors` (an object equal to theirs; dev/plan.md 2.6b, G5). False if not
  // possible (the default).
  virtual bool repoint_to(const Factors* factors) { return false; }
  
  virtual void close_ofstreams()=0;
  
  virtual Parameters get_current_algorithm_parameters() const=0;
  
protected:
  
  void make_copy(const MoveOutput &another);

};
}

#endif
