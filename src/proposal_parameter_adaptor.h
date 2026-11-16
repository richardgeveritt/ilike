#ifndef PROPOSALPARAMETERADAPTOR_H
#define PROPOSALPARAMETERADAPTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <string>

#include "mcmc_adaptor.h"
#include "ilike_header.h"

namespace ilike
{

// The general controlled MCMC of Andrieu & Thoms (2008), Algorithm 1, in the
// stochastic approximation form (23): theta_{i+1} = theta_i + gamma_{i+1} H,
// for named proposal parameters theta of a proposal written by the user
// (proposal_parameters in its functions; dev/plan.md 3.2). Two forms:
//   ACCEPTANCE_RATE  for each named parameter (positive), on the log scale:
//                    log theta += gamma (alpha - target)  (eq. (20));
//   CUSTOM           H returned by the user's function
//                    adapt_m_proposal / adapt_mh_proposal
//                    (proposal_parameters, parameters = X_i,
//                     proposed_parameters = Y_{i+1},
//                     new_parameters = X_{i+1}, acceptance probability):
//                    every parameter named in H is updated.
// Each chain adapts its own copy of the proposal parameters: the adaptor
// holds them, and points its kernel at them before every move.
// With averaging (Polyak & Juditsky 1992) the recursion runs as above and
// the kernel uses the average of its iterates (for ACCEPTANCE_RATE, of the
// log parameters); a user's H is given the recursion's own iterate.
class ProposalParameterAdaptor : public MCMCAdaptor
{

public:

  ProposalParameterAdaptor();

  // acceptance-rate form
  ProposalParameterAdaptor(const std::vector<std::string> &parameter_names_in,
                           double target_acceptance_in,
                           double gain_exponent_in,
                           double gain_constant_in,
                           bool averaging_in = false);

  // custom form
  ProposalParameterAdaptor(AdaptProposalPtr adapt_in,
                           double gain_exponent_in,
                           double gain_constant_in,
                           bool averaging_in = false);

  virtual ~ProposalParameterAdaptor();

  ProposalParameterAdaptor(const ProposalParameterAdaptor &another);
  void operator=(const ProposalParameterAdaptor &another);
  MCMCAdaptor* duplicate() const;

  void prepare_for_move();
  void adapt(const AdaptationStep &step);
  void append_state(std::vector<std::string> *names,
                    std::vector<double> *values,
                    bool final,
                    const std::string &prefix) const;

protected:

  void make_copy(const ProposalParameterAdaptor &another);

  AdaptProposalPtr adapt_function;
  std::vector<std::string> parameter_names; // acceptance-rate form
  double target_acceptance;
  double gain_exponent;
  double gain_constant;
  size_t number_of_updates;
  Parameters own_parameters; // this chain's proposal parameters (what the kernel uses)
  Parameters iterate; // the recursion's theta_n (log scale for ACCEPTANCE_RATE)
  Parameters average; // its running average (with averaging)
  bool averaging;
  bool started;

};

}

#endif
