#ifndef ELEMENTRANDOMWALKPROPOSALKERNEL_H
#define ELEMENTRANDOMWALKPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <string>
#include "symmetric_proposal_kernel.h"
#include "particle.h"
#include "distributions.h"
#include "ilike_header.h"

namespace ilike
{

// how an element moves: a normal step, a uniform step, or 1 - x (a 0/1 element)
// normal_conditional: no step, a draw from the element's normal full
// conditional (gibbs on a normal mean, C206; only in InPlaceElementScanMCMC);
// gamma_conditional, invgamma_conditional, beta_conditional: likewise from a
// conjugate gamma, inverse gamma or beta full conditional (C215)
enum class ElementStep { normal, uniform, flip, normal_conditional, gamma_conditional, invgamma_conditional, beta_conditional };

// A random walk on one element of a variable (dev/plan.md I4c): the element
// (row, column, 0-based) moves by a normal draw with standard deviation
// `scale` (normal), a uniform draw on [-scale, scale] (uniform), or becomes
// 1 - x (flip, for a 0/1 element; symmetric); the other elements stay. Written for scans over the elements of a variable
// (`for (j in 1:J) move m = metropolis(norm_rw(theta[j], s))`).
class ElementRandomWalkProposalKernel : public SymmetricProposalKernel
{
public:

  ElementRandomWalkProposalKernel();
  ElementRandomWalkProposalKernel(const std::string &variable_in,
                                  size_t row_in,
                                  size_t column_in,
                                  double scale_in,
                                  ElementStep step_in);
  virtual ~ElementRandomWalkProposalKernel();

  ElementRandomWalkProposalKernel(const ElementRandomWalkProposalKernel &another);
  void operator=(const ElementRandomWalkProposalKernel &another);
  Kernel* duplicate() const;
  ProposalKernel* proposal_kernel_duplicate() const;
  SymmetricProposalKernel* symmetric_proposal_kernel_duplicate() const;

  void set_proposal_parameters(Parameters* proposal_parameters_in);
  GradientEstimatorOutput* simulate_gradient_estimator_output() const;
  std::vector<const ProposalKernel*> get_proposals() const;
  void set_index(Index* index_in);
  void set_index_if_null(Index* index_in);
  bool can_be_evaluated() const;
  void set_data(Data* data_in);

  // the step (a normal sd or a uniform half-width), adapted by an
  // ElementScaleAdaptor or an ElementSMCAdaptor; and the element moved
  double get_scale() const { return this->scale; }
  void set_scale(double scale_in) { this->scale = scale_in; }
  const std::string &get_variable() const { return this->variable; }
  size_t get_row() const { return this->row; }
  size_t get_column() const { return this->column; }
  ElementStep get_step_kind() const { return this->step_kind; }

protected:

  double specific_evaluate_kernel(const Particle &proposed_particle,
                                  const Particle &old_particle) const;
  double specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                            const Particle &old_particle) const;
  Parameters simulate(RandomNumberGenerator &rng,
                      const Particle &particle) const;
  void simulate_into(RandomNumberGenerator &rng,
                     const Particle &particle,
                     ProposalTarget &target) const;
  Parameters simulate(RandomNumberGenerator &rng,
                      const Particle &particle,
                      const Parameters &conditioned_on_parameters) const;
  Parameters subsample_simulate(RandomNumberGenerator &rng,
                                const Particle &particle) const;
  Parameters subsample_simulate(RandomNumberGenerator &rng,
                                const Particle &particle,
                                const Parameters &conditioned_on_parameters) const;
  Parameters subsample_simulate(RandomNumberGenerator &rng,
                                const std::string &variable,
                                const Particle &particle) const;
  Parameters subsample_simulate(RandomNumberGenerator &rng,
                                const std::string &variable,
                                const Particle &particle,
                                const Parameters &conditioned_on_parameters) const;
  arma::mat specific_gradient_of_log(const std::string &variable,
                                     const Particle &proposed_particle,
                                     const Particle &old_particle);
  arma::mat specific_subsample_gradient_of_log(const std::string &variable,
                                               const Particle &proposed_particle,
                                               const Particle &old_particle);

  // the element's new value given its current one
  double step(RandomNumberGenerator &rng, double current) const;
  const arma::mat &current_value(const Particle &particle) const;

  void make_copy(const ElementRandomWalkProposalKernel &another);

  std::string variable;
  size_t row;
  size_t column;
  double scale;
  ElementStep step_kind;
};

}

#endif
