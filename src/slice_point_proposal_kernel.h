#ifndef SLICEPOINTPROPOSALKERNEL_H
#define SLICEPOINTPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <string>
#include <vector>
#include "proposal_kernel.h"
#include "particle.h"
#include "distributions.h"
#include "ilike_header.h"

namespace ilike
{

// The settings of a slice sampling move, move m = slice(x, width, lower,
// upper, max_steps) (dev/plan.md 3.21 (c3); Neal 2003), and the means of
// building its candidate points: a particle with element `element` of x set
// to `value` (SliceMCMC sets them, then builds the particle through
// ProposalKernel::move, which gives it new outputs for the move's factors).
// Not a proposal with a density: the recipe makes the move a SliceMCMC.
class SlicePointProposalKernel : public ProposalKernel
{
public:

  SlicePointProposalKernel();
  virtual ~SlicePointProposalKernel();
  SlicePointProposalKernel(const std::string &variable_in, double width_in,
                           double lower_in, double upper_in, size_t max_steps_in);

  const std::string &get_variable() const { return this->variable; }
  double get_width() const { return this->width; }
  double get_lower() const { return this->lower; }
  double get_upper() const { return this->upper; }
  size_t get_max_steps() const { return this->max_steps; }
  // the candidate's element and value (read by simulate)
  void set_point(arma::uword element_in, double value_in) const { this->element = element_in; this->value = value_in; }
  SlicePointProposalKernel(const SlicePointProposalKernel &another);
  void operator=(const SlicePointProposalKernel &another);
  Kernel* duplicate() const;
  ProposalKernel* proposal_kernel_duplicate() const;

  void set_proposal_parameters(Parameters* proposal_parameters_in);
  GradientEstimatorOutput* simulate_gradient_estimator_output() const;
  std::vector<const ProposalKernel*> get_proposals() const;
  void set_index(Index* index_in);
  void set_index_if_null(Index* index_in);
  bool can_be_evaluated() const;
  void set_data(Data* data_in);

protected:

  double specific_evaluate_kernel(const Particle &proposed_particle,
                                  const Particle &old_particle) const;
  double specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                            const Particle &old_particle) const;
  Parameters simulate(RandomNumberGenerator &rng,
                      const Particle &particle) const;
  Parameters subsample_simulate(RandomNumberGenerator &rng,
                                const Particle &particle) const;
  Parameters subsample_simulate(RandomNumberGenerator &rng,
                                const std::string &variable,
                                const Particle &particle) const;
  arma::mat specific_gradient_of_log(const std::string &variable,
                                     const Particle &proposed_particle,
                                     const Particle &old_particle);
  arma::mat specific_subsample_gradient_of_log(const std::string &variable,
                                               const Particle &proposed_particle,
                                               const Particle &old_particle);

  void make_copy(const SlicePointProposalKernel &another);
  std::string variable;
  double width;
  double lower;
  double upper;
  size_t max_steps;
  mutable arma::uword element = 0;
  mutable double value = 0.0;
};

}

#endif
