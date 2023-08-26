#ifndef UNIFORMRANDOMWALKPROPOSALKERNEL_H
#define UNIFORMRANDOMWALKPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <boost/unordered_map.hpp>

#include "symmetric_proposal_kernel.h"
#include "particle.h"
#include "distributions.h"
#include "ilike_header.h"

class UniformRandomWalkProposalKernel : public SymmetricProposalKernel
{

public:

  UniformRandomWalkProposalKernel();
  virtual ~UniformRandomWalkProposalKernel();
  
  // make cov_names from var_names and find cov adaptively
  UniformRandomWalkProposalKernel(const std::vector<std::string> &variable_names_in);
  
  // make cov_names from var_names
  UniformRandomWalkProposalKernel(const std::vector<std::string> &variable_names_in,
                                  const std::vector<arma::mat> &halfwidths_in);
  
  UniformRandomWalkProposalKernel(const std::string &variable_name_in,
                                  const double &halfwidth_in);

  UniformRandomWalkProposalKernel(const std::string &variable_name_in,
                                  const arma::colvec &halfwidth_in);
  
  UniformRandomWalkProposalKernel(const std::string &variable_name_in,
                                  const arma::mat &halfwidth_in);

  UniformRandomWalkProposalKernel(const UniformRandomWalkProposalKernel &another);

  void operator=(const UniformRandomWalkProposalKernel &another);
  Kernel* duplicate() const;
  ProposalKernel* proposal_kernel_duplicate() const;
  SymmetricProposalKernel* symmetric_proposal_kernel_duplicate() const;
  
  void set_proposal_parameters(Parameters* proposal_parameters_in);
  
  GradientEstimatorOutput* simulate_gradient_estimator_output() const;
  
  std::vector<const ProposalKernel*> get_proposals() const;
  
  void set_index(Index* index_in);
  
// Mh has its own parameters.
  // Stochastic has some weights.
  // MH has sim prop and eval prop, take in params. Use current value in acceptance, Set current value if accepted.
  // Proposal needs to call simulate in all llhdoutputs

protected:
  
  double specific_evaluate_kernel(const Particle &proposed_particle,
                                  const Particle &old_particle) const;
  /*
  double specific_evaluate_kernel(Particle &proposed_particle,
                                  Particle &old_particle,
                                  const Parameters &conditioned_on_parameters) const;
  */
  double specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                            const Particle &old_particle) const;
  
  /*
  double specific_subsample_evaluate_kernel(Particle &proposed_particle,
                                            Particle &old_particle,
                                            const Parameters &conditioned_on_parameters) const;
  */
  
  Parameters simulate(RandomNumberGenerator &rng,
                      const Particle &particle) const;
  
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
  
  /*
  arma::mat specific_gradient_of_log(const std::string &variable,
                                     Particle &proposed_particle,
                                     Particle &old_particle,
                                     const Parameters &conditioned_on_parameters);
  */
  
  arma::mat specific_subsample_gradient_of_log(const std::string &variable,
                                               const Particle &proposed_particle,
                                               const Particle &old_particle);
  
  /*
  arma::mat specific_subsample_gradient_of_log(const std::string &variable,
                                               Particle &proposed_particle,
                                               Particle &old_particle,
                                               const Parameters &conditioned_on_parameters);
  */

  void make_copy(const UniformRandomWalkProposalKernel &another);
  
  //bool unused_variables_kept;
  
  // store the half widths for every variable
  boost::unordered_map< std::string, arma::mat> proposal_info;

};

#endif
