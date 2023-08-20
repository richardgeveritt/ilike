#ifndef PDMPMCMCOUTPUT_H
#define PDMPMCMCOUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "move_output.h"

class PDMPMCMCOutput : public MoveOutput
{

public:

  PDMPMCMCOutput();

  virtual ~PDMPMCMCOutput();

  PDMPMCMCOutput(const PDMPMCMCOutput &another);

  void operator=(const PDMPMCMCOutput &another);
  MoveOutput* duplicate() const;

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
  
  Particle dummy;

  void make_copy(const PDMPMCMCOutput &another);

};

#endif
