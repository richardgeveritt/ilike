#ifndef GAUSSIANINDEPENDENTPROPOSALKERNEL_H
#define GAUSSIANINDEPENDENTPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "independent_proposal_kernel.h"
#include "particle.h"
#include "distributions.h"
#include "ilike_header.h"
#include "gaussian_proposal_info.h"

class GaussianIndependentProposalKernel : public IndependentProposalKernel
{

public:

  GaussianIndependentProposalKernel();
  virtual ~GaussianIndependentProposalKernel();
  
  // find mean and cov adaptively
  GaussianIndependentProposalKernel(const std::vector<std::string> &variable_names_in);
  
  // find cov adaptively
  GaussianIndependentProposalKernel(const std::vector<std::string> &variable_names_in,
                                    const std::vector<arma::colvec> &means_in);
  
  // find mean adaptively
  GaussianIndependentProposalKernel(const std::vector<std::string> &variable_names_in,
                                    const std::vector<arma::mat> &covariances_in);
  
  GaussianIndependentProposalKernel(const std::string &variable_name_in,
                                    const double &mean_in,
                                    const double &sd_in);
  
  GaussianIndependentProposalKernel(const std::vector<std::string> &variable_names_in,
                                    const std::vector<arma::colvec> &means_in,
                                    const std::vector<arma::mat> &covariances_in);

  GaussianIndependentProposalKernel(const GaussianIndependentProposalKernel &another);

  void operator=(const GaussianIndependentProposalKernel &another);
  Kernel* duplicate() const;
  ProposalKernel* proposal_kernel_duplicate() const;
  IndependentProposalKernel* independent_proposal_kernel_duplicate() const;
  
  void set_mean(const std::string &variable,
                const arma::colvec &mean_in);
  
  void set_covariance(const std::string &variable,
                      const arma::mat &covariance_in);
  
  arma::mat get_inverse_covariance(const std::string &variable);
  arma::mat get_covariance(const std::string &variable);
  arma::mat get_covariance(const std::vector<std::string> &variables);
  
  double evaluate_independent_kernel(const Parameters &proposed_particle) const;
  //double evaluate_independent_kernel(Variables* proposed_particle,
  //                                   const Parameters &conditioned_on_parameters) const;
  double subsample_evaluate_independent_kernel(const Parameters &proposed_particle) const;
  
  Parameters independent_simulate(RandomNumberGenerator &rng) const;
  
  Parameters independent_simulate(RandomNumberGenerator &rng,
                                  const Parameters &conditioned_on_parameters) const;
  
  Parameters subsample_independent_simulate(RandomNumberGenerator &rng) const;
  
  Parameters subsample_independent_simulate(RandomNumberGenerator &rng,
                                            const Parameters &conditioned_on_parameters) const;
  
  Parameters subsample_independent_simulate(RandomNumberGenerator &rng,
                                            const std::string &variable) const;
  
  Parameters subsample_independent_simulate(RandomNumberGenerator &rng,
                                            const std::string &variable,
                                            const Parameters &conditioned_on_parameters) const;
  
  arma::mat independent_gradient_of_log(const std::string &variable,
                                        const Parameters &proposed_particle);
  arma::mat subsample_independent_gradient_of_log(const std::string &variable,
                                                  const Parameters &proposed_particle);
  
// Mh has its own parameters.
  // Stochastic has some weights.
  // MH has sim prop and eval prop, t**ake in params. Use current value in acceptance, Set current value if accepted.
  // Proposal needs to call simulate in all llhdoutputs

protected:

  void make_copy(const GaussianIndependentProposalKernel &another);
  
  boost::unordered_map< std::string, GaussianProposalInfo> proposal_info;
  
};

#endif
