#ifndef MOVEOUTPUT_H
#define MOVEOUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <iostream>
#include "particle.h"
#include "ilike_header.h"

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
  
  virtual void write_factors(const std::string &directory_name,
                             const std::string &index) const=0;
  virtual void write_ensemble_factors(const std::string &directory_name,
                                      const std::string &index) const=0;
  
  virtual size_t length() const=0;
  
  virtual void close_ofstreams()=0;

protected:

  void make_copy(const MoveOutput &another);
  
  //std::vector<std::string> vector_variables;
  //std::vector<std::string> any_variables;

};

#endif
