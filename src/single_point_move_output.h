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

  SinglePointMoveOutput(const Parameters &parameter_in,
                        const Factors* factors_in,
                        const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
                        const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in);
  SinglePointMoveOutput(const Parameters &parameter_in,
                        const EnsembleFactors* factors_in);
  
  SinglePointMoveOutput(Parameters &&parameter_in,
                        const Factors* factors_in,
                        const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
                        const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in);
  SinglePointMoveOutput(Parameters &&parameter_in,
                        const EnsembleFactors* factors_in);
  
  SinglePointMoveOutput(const Particle &particle_in);
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
  
  size_t length() const;
  
  void close_ofstreams();
  
  Parameters get_current_algorithm_parameters() const;
  
protected:

  void make_copy(const SinglePointMoveOutput &another);
  
  Particle output;
  
  Parameters algorithm_parameters;

};

#endif
