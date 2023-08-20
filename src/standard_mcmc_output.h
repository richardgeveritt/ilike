#ifndef STANDARDMCMCOUTPUT_H
#define STANDARDMCMCOUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "move_output.h"
#include "parameters.h"
#include "ilike_header.h"

class StandardMCMCOutput : public MoveOutput
{

public:

  //StandardMCMCOutput(const Parameters &parameter_in);
  StandardMCMCOutput();

  virtual ~StandardMCMCOutput();

  StandardMCMCOutput(const StandardMCMCOutput &another);

  void operator=(const StandardMCMCOutput &another);
  MoveOutput* duplicate() const;
  
  void push_back(const Particle &particle_in);

  Particle& back();
  Particle back() const;
  
  std::vector<Parameters> get_vector_of_parameters() const;
  
  void write_vector_points(const std::vector<std::string> &variables,
                           std::ofstream &file_stream,
                           std::shared_ptr<Transform> transform) const;
  void write_any_points(const std::vector<std::string> &variables,
                        std::ofstream &file_stream) const;
  
  void write_factors(const std::string &directory_name,
                     const std::string &index) const;
  void write_ensemble_factors(const std::string &directory_name,
                              const std::string &index) const;
  
  size_t length() const;
  
  void close_ofstreams();
  
protected:

  void make_copy(const StandardMCMCOutput &another);
  
  std::deque<Particle> output;

};

#endif
