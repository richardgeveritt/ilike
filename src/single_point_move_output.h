#ifndef SINGLEPOINTMOVEOUTPUT_H
#define SINGLEPOINTMOVEOUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "move_output.h"
#include "parameters.h"
#include "particle.h"
#include "ilike_header.h"

class Factors;

class SinglePointMoveOutput : public MoveOutput
{

public:

  SinglePointMoveOutput(Parameters &&parameter_in,
                        Factors* factors_in);
  SinglePointMoveOutput(Parameters &&parameter_in,
                        EnsembleFactors* factors_in);
  SinglePointMoveOutput(const Parameters &parameter_in,
                        EnsembleFactors* factors_in);
  //SinglePointMoveOutput(const Particle &particle_in);
  SinglePointMoveOutput(Particle &&particle_in);
  SinglePointMoveOutput();

  virtual ~SinglePointMoveOutput();

  SinglePointMoveOutput(const SinglePointMoveOutput &another);

  void operator=(const SinglePointMoveOutput &another);
  MoveOutput* duplicate() const;

  Particle& back();
  Particle back() const;
  
  std::vector<Parameters> get_vector_of_parameters() const;
  
  void write_vector_points(const std::vector<std::string> &variables,
                           std::ofstream &file_stream,
                           std::shared_ptr<Transform> inverse_transform) const;
  void write_any_points(const std::vector<std::string> &variables,
                        std::ofstream &file_stream) const;
  
  void write_factors(const std::string &directory_name,
                     const std::string &index) const;
  void write_ensemble_factors(const std::string &directory_name,
                              const std::string &index) const;
  
  void close_ofstreams();
  
protected:

  void make_copy(const SinglePointMoveOutput &another);
  
  Particle output;

};

#endif
