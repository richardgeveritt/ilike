#include <cmath>
#include <algorithm>
#include "proposal_parameter_adaptor.h"
#include "proposal_kernel.h"
#include "particle.h"

namespace ilike
{

static bool has_parameter(const Parameters &parameters, const std::string &name)
{
  std::vector<std::string> names = parameters.get_vector_variables();
  return std::find(names.begin(), names.end(), name)!=names.end();
}

ProposalParameterAdaptor::ProposalParameterAdaptor()
:MCMCAdaptor()
{
  this->adapt_function = NULL;
  this->target_acceptance = 0.234;
  this->gain_exponent = 0.6;
  this->gain_constant = 1.0;
  this->number_of_updates = 0;
  this->averaging = false;
  this->started = false;
}

ProposalParameterAdaptor::ProposalParameterAdaptor(const std::vector<std::string> &parameter_names_in,
                                                   double target_acceptance_in,
                                                   double gain_exponent_in,
                                                   double gain_constant_in,
                                                   bool averaging_in)
:MCMCAdaptor()
{
  this->averaging = averaging_in;
  if ( (target_acceptance_in<=0.0) || (target_acceptance_in>=1.0) )
    Rcpp::stop("MCMC adaptor: the target acceptance rate must be between 0 and 1.");
  if ( (gain_exponent_in<=0.5) || (gain_exponent_in>1.0) )
    Rcpp::stop("MCMC adaptor: the gain exponent must be in (0.5, 1].");
  this->adapt_function = NULL;
  this->parameter_names = parameter_names_in;
  this->target_acceptance = target_acceptance_in;
  this->gain_exponent = gain_exponent_in;
  this->gain_constant = gain_constant_in;
  this->number_of_updates = 0;
  this->started = false;
}

ProposalParameterAdaptor::ProposalParameterAdaptor(AdaptProposalPtr adapt_in,
                                                   double gain_exponent_in,
                                                   double gain_constant_in,
                                                   bool averaging_in)
:MCMCAdaptor()
{
  this->averaging = averaging_in;
  if ( (gain_exponent_in<=0.5) || (gain_exponent_in>1.0) )
    Rcpp::stop("MCMC adaptor: the gain exponent must be in (0.5, 1].");
  this->adapt_function = adapt_in;
  this->target_acceptance = 0.0;
  this->gain_exponent = gain_exponent_in;
  this->gain_constant = gain_constant_in;
  this->number_of_updates = 0;
  this->started = false;
}

ProposalParameterAdaptor::~ProposalParameterAdaptor()
{
}

ProposalParameterAdaptor::ProposalParameterAdaptor(const ProposalParameterAdaptor &another)
:MCMCAdaptor(another)
{
  this->make_copy(another);
}

void ProposalParameterAdaptor::operator=(const ProposalParameterAdaptor &another)
{
  if(this == &another)
    return;
  MCMCAdaptor::operator=(another);
  this->make_copy(another);
}

void ProposalParameterAdaptor::make_copy(const ProposalParameterAdaptor &another)
{
  this->adapt_function = another.adapt_function;
  this->parameter_names = another.parameter_names;
  this->target_acceptance = another.target_acceptance;
  this->gain_exponent = another.gain_exponent;
  this->gain_constant = another.gain_constant;
  this->number_of_updates = another.number_of_updates;
  // deep copies: Parameters share the matrices of "fixed" entries, and each
  // chain's adaptor must have its own
  this->own_parameters = another.own_parameters.deep_copy();
  this->iterate = another.iterate.deep_copy();
  this->average = another.average.deep_copy();
  this->averaging = another.averaging;
  this->started = another.started;
}

MCMCAdaptor* ProposalParameterAdaptor::duplicate() const
{
  return( new ProposalParameterAdaptor(*this));
}

void ProposalParameterAdaptor::prepare_for_move()
{
  if (this->proposal==NULL)
    Rcpp::stop("MCMC adaptor: not attached to a proposal.");
  if (!this->started)
  {
    // theta_0: the proposal parameters the algorithm was given
    Parameters* given = this->proposal->get_proposal_parameters();
    if (given==NULL)
      Rcpp::stop("MCMC adaptor: adapting named proposal parameters needs a proposal written by the user that takes proposal_parameters.");
    // a deep copy: the algorithm's parameters (shared by every chain) mark
    // their entries as fixed, and a plain copy would share their matrices
    this->own_parameters = given->deep_copy();
    for (auto i=this->parameter_names.begin(); i!=this->parameter_names.end(); ++i)
    {
      if (!has_parameter(this->own_parameters, *i))
        Rcpp::stop("MCMC adaptor: proposal parameter " + *i + " is not given (algorithm_parameter_list).");
      if (arma::any(arma::vectorise(this->own_parameters[*i])<=0.0))
        Rcpp::stop("MCMC adaptor: ilike::acceptance_rate adapts positive parameters (on the log scale); " + *i + " is not positive.");
    }
    // theta_0, on the scale of the recursion
    this->iterate = this->own_parameters.deep_copy();
    if (this->adapt_function==NULL)
    {
      for (auto i=this->parameter_names.begin(); i!=this->parameter_names.end(); ++i)
        this->iterate[*i] = arma::log(this->own_parameters[*i]);
    }
    this->average = this->iterate.deep_copy();
    this->started = true;
  }
  // this chain's own copy (a copied kernel still points at the parameters
  // of the kernel it was copied from)
  this->proposal->set_proposal_parameters(&this->own_parameters);
}

void ProposalParameterAdaptor::adapt(const AdaptationStep &step)
{
  this->number_of_updates = this->number_of_updates + 1;
  double gain = this->gain_constant*pow(double(this->number_of_updates+1), -this->gain_exponent);
  std::vector<std::string> updated;
  if (this->adapt_function==NULL)
  {
    // log theta += gamma (alpha - target)
    for (auto i=this->parameter_names.begin(); i!=this->parameter_names.end(); ++i)
      this->iterate[*i] = this->iterate[*i] + gain*(step.acceptance_probability - this->target_acceptance);
    updated = this->parameter_names;
  }
  else
  {
    Parameters H = this->adapt_function(this->iterate,
                                        step.previous->parameters,
                                        step.proposed->parameters,
                                        step.current->parameters,
                                        step.acceptance_probability);
    std::vector<std::string> names = H.get_vector_variables();
    for (auto i=names.begin(); i!=names.end(); ++i)
    {
      if (!has_parameter(this->own_parameters, *i))
        Rcpp::stop("MCMC adaptor: the adaptation function returned " + *i + ", which is not a proposal parameter.");
      arma::mat current = this->iterate[*i];
      arma::mat increment = H[*i];
      if ( (current.n_rows!=increment.n_rows) || (current.n_cols!=increment.n_cols) )
        Rcpp::stop("MCMC adaptor: the adaptation function returned " + *i + " with the wrong dimensions.");
      this->iterate[*i] = current + gain*increment;
    }
    updated = this->own_parameters.get_vector_variables();
  }
  // theta-bar_n = theta-bar_{n-1} + (theta_n - theta-bar_{n-1})/(n+1)
  double weight = 1.0/double(this->number_of_updates + 1);
  for (auto i=updated.begin(); i!=updated.end(); ++i)
  {
    arma::mat used;
    if (this->averaging)
    {
      this->average[*i] = this->average[*i] + weight*(this->iterate[*i] - this->average[*i]);
      used = this->average[*i];
    }
    else
    {
      used = this->iterate[*i];
    }
    if (this->adapt_function==NULL)
      used = arma::exp(used);
    this->own_parameters[*i] = used;
  }
}

void ProposalParameterAdaptor::append_state(std::vector<std::string> *names,
                                            std::vector<double> *values,
                                            bool final,
                                            const std::string &prefix) const
{
  if (!this->started)
    return;
  // the parameters adapted by name; with an adaptation function, all the
  // proposal parameters (so the columns do not depend on what H returns)
  std::vector<std::string> written = this->parameter_names;
  if (this->adapt_function!=NULL)
    written = this->own_parameters.get_vector_variables();
  for (auto i=written.begin(); i!=written.end(); ++i)
  {
    arma::colvec v = arma::vectorise(this->own_parameters[*i]);
    for (size_t k=0; k<v.n_elem; ++k)
    {
      if (names!=NULL)
        names->push_back(prefix + *i + "[" + std::to_string(k+1) + "]");
      if (values!=NULL)
        values->push_back(v[k]);
    }
  }
}

}
