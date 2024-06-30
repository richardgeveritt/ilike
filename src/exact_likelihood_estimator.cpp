#include "exact_likelihood_estimator.h"
#include "exact_likelihood_estimator_output.h"
#include "independent_proposal_kernel.h"
#include "distribution_factor.h"
#include "likelihood_factor.h"
#include "custom_distribution_factor.h"
#include "custom_likelihood_factor.h"

ExactLikelihoodEstimator::ExactLikelihoodEstimator()
:LikelihoodEstimator()
{
}

ExactLikelihoodEstimator::ExactLikelihoodEstimator(RandomNumberGenerator* rng_in,
                                                   size_t* seed_in,
                                                   Data* data_in,
                                                   EvaluateLogLikelihoodPtr llhd_in,
                                                   bool smcfixed_flag_in)
:LikelihoodEstimator(rng_in, seed_in, data_in, Parameters(), smcfixed_flag_in)
{
  this->numerator_likelihood_factors.push_back(new CustomLikelihoodFactor(llhd_in,
                                                                          data_in));
  //this->output = new ExactLikelihoodEstimatorOutput();
}

ExactLikelihoodEstimator::ExactLikelihoodEstimator(RandomNumberGenerator* rng_in,
                                                   size_t* seed_in,
                                                   Data* data_in,
                                                   EvaluateLogDistributionPtr dist_in,
                                                   bool smcfixed_flag_in)
:LikelihoodEstimator(rng_in, seed_in, data_in, Parameters(), smcfixed_flag_in)
{
  this->numerator_distribution_factors.push_back(new CustomDistributionFactor(dist_in));
  //this->output = new ExactDistributionEstimatorOutput();
}

ExactLikelihoodEstimator::ExactLikelihoodEstimator(RandomNumberGenerator* rng_in,
                                                   size_t* seed_in,
                                                   Data* data_in,
                                                   DistributionFactor* dist_in,
                                                   bool smcfixed_flag_in)
:LikelihoodEstimator(rng_in, seed_in, data_in, Parameters(), smcfixed_flag_in)
{
  this->numerator_distribution_factors.push_back(dist_in);
}

ExactLikelihoodEstimator::ExactLikelihoodEstimator(RandomNumberGenerator* rng_in,
                                                   size_t* seed_in,
                                                   Data* data_in,
                                                   const std::vector<DistributionFactor*> &numerator_distribution_factors_in,
                                                   bool smcfixed_flag_in)
:LikelihoodEstimator(rng_in, seed_in, data_in, Parameters(), smcfixed_flag_in)
{
  this->numerator_distribution_factors = numerator_distribution_factors_in;
}

ExactLikelihoodEstimator::ExactLikelihoodEstimator(RandomNumberGenerator* rng_in,
                                                   size_t* seed_in,
                                                   Data* data_in,
                                                   LikelihoodFactor* llhd_in,
                                                   bool smcfixed_flag_in)
:LikelihoodEstimator(rng_in, seed_in, data_in, Parameters(), smcfixed_flag_in)
{
  this->numerator_likelihood_factors.push_back(llhd_in);
}

ExactLikelihoodEstimator::ExactLikelihoodEstimator(RandomNumberGenerator* rng_in,
                                                   size_t* seed_in,
                                                   Data* data_in,
                                                   const std::vector<LikelihoodFactor*> &numerator_likelihood_factors_in,
                                                   bool smcfixed_flag_in)
:LikelihoodEstimator(rng_in, seed_in, data_in, Parameters(), smcfixed_flag_in)
{
  this->numerator_likelihood_factors = numerator_likelihood_factors_in;
}

ExactLikelihoodEstimator::ExactLikelihoodEstimator(RandomNumberGenerator* rng_in,
                                                   size_t* seed_in,
                                                   Data* data_in,
                                                   const std::vector<DistributionFactor*> &numerator_distribution_factors_in,
                                                   const std::vector<LikelihoodFactor*> &numerator_likelihood_factors_in,
                                                   bool smcfixed_flag_in)
:LikelihoodEstimator(rng_in, seed_in, data_in, Parameters(), smcfixed_flag_in)
{
  this->numerator_distribution_factors = numerator_distribution_factors_in;
  this->numerator_likelihood_factors = numerator_likelihood_factors_in;
}

ExactLikelihoodEstimator::ExactLikelihoodEstimator(RandomNumberGenerator* rng_in,
                                                   size_t* seed_in,
                                                   Data* data_in,
                                                   const std::vector<DistributionFactor*> &numerator_distribution_factors_in,
                                                   const std::vector<LikelihoodFactor*> &numerator_likelihood_factors_in,
                                                   const std::vector<ProposalKernel*> &numerator_likelihood_proposals_in,
                                                   bool smcfixed_flag_in)
:LikelihoodEstimator(rng_in, seed_in, data_in, Parameters(), smcfixed_flag_in)
{
  this->numerator_distribution_factors = numerator_distribution_factors_in;
  this->numerator_likelihood_factors = numerator_likelihood_factors_in;
  this->numerator_likelihood_proposals = numerator_likelihood_proposals_in;
}

ExactLikelihoodEstimator::ExactLikelihoodEstimator(RandomNumberGenerator* rng_in,
                                                   size_t* seed_in,
                                                   Data* data_in,
                                                   IndependentProposalKernel* dist_in,
                                                   bool smcfixed_flag_in)
:LikelihoodEstimator(rng_in, seed_in, data_in, Parameters(), smcfixed_flag_in)
{
  this->numerator_proposals.push_back(dist_in);
}

ExactLikelihoodEstimator::ExactLikelihoodEstimator(RandomNumberGenerator* rng_in,
                                                   size_t* seed_in,
                                                   Data* data_in,
                                                   ProposalKernel* dist_in,
                                                   bool smcfixed_flag_in)
:LikelihoodEstimator(rng_in, seed_in, data_in, Parameters(), smcfixed_flag_in)
{
  this->numerator_likelihood_proposals.push_back(dist_in);
}

ExactLikelihoodEstimator::ExactLikelihoodEstimator(RandomNumberGenerator* rng_in,
                                                   size_t* seed_in,
                                                   Data* data_in,
                                                   EvaluateLogDistributionPtr prior_in,
                                                   EvaluateLogLikelihoodPtr llhd_in,
                                                   bool smcfixed_flag_in)
:LikelihoodEstimator(rng_in, seed_in, data_in, Parameters(), smcfixed_flag_in)
{
  this->numerator_likelihood_factors.push_back(new CustomLikelihoodFactor(llhd_in,
                                                                          data_in));
  
  this->numerator_distribution_factors.push_back(new CustomDistributionFactor(prior_in));
}

/*
 ExactLikelihoodEstimator::ExactLikelihoodEstimator(RandomNumberGenerator* rng_in,
 size_t* seed_in,
 const Data* data_in,
 EvaluateLogDistributionPtr prior_in,
 EvaluateLogLikelihoodPtr llhd_in,
 EvaluateLogDistributionPtr proposal_in,
 bool smcfixed_flag_in)
 :LikelihoodEstimator(rng_in, seed_in, data_in)
 {
 this->numerator_llhds.reserve(1);
 this->numerator_llhds.push_back(llhd_in);
 
 this->numerator_distributions.reserve(1);
 this->numerator_distributions.push_back(prior_in);
 
 this->denominator_distributions.reserve(1);
 this->denominator_distributions.push_back(proposal_in);
 
 this->smcfixed_flag = smcfixed_flag_in;
 }
 */

ExactLikelihoodEstimator::~ExactLikelihoodEstimator()
{
  for (auto i=this->numerator_distribution_factors.begin();
       i!=this->numerator_distribution_factors.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  
  for (auto i=this->numerator_likelihood_factors.begin();
       i!=this->numerator_likelihood_factors.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  
  /*
  for (auto i=this->numerator_proposals.begin();
       i!=this->numerator_proposals.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  */
}

//Copy constructor for the ExactLikelihoodEstimator class.
ExactLikelihoodEstimator::ExactLikelihoodEstimator(const ExactLikelihoodEstimator &another)
:LikelihoodEstimator(another)
{
  this->make_copy(another);
}

void ExactLikelihoodEstimator::operator=(const ExactLikelihoodEstimator &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  for (auto i=this->numerator_distribution_factors.begin();
       i!=this->numerator_distribution_factors.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  this->numerator_distribution_factors.clear();
  
  for (auto i=this->numerator_likelihood_factors.begin();
       i!=this->numerator_likelihood_factors.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  this->numerator_likelihood_factors.clear();
  
  LikelihoodEstimator::operator=(another);
  this->make_copy(another);
}

LikelihoodEstimator* ExactLikelihoodEstimator::duplicate() const
{
  return( new ExactLikelihoodEstimator(*this));
}

void ExactLikelihoodEstimator::make_copy(const ExactLikelihoodEstimator &another)
{
  this->numerator_distribution_factors.resize(0);
  this->numerator_distribution_factors.reserve(another.numerator_distribution_factors.size());
  for (auto i=another.numerator_distribution_factors.begin();
       i!=numerator_distribution_factors.end();
       ++i)
  {
    if (*i!=NULL)
      this->numerator_distribution_factors.push_back((*i)->distribution_factor_duplicate());
    else
      this->numerator_distribution_factors.push_back(NULL);
  }
  
  this->numerator_likelihood_factors.resize(0);
  this->numerator_likelihood_factors.reserve(another.numerator_likelihood_factors.size());
  for (auto i=another.numerator_likelihood_factors.begin();
       i!=numerator_likelihood_factors.end();
       ++i)
  {
    if (*i!=NULL)
      this->numerator_likelihood_factors.push_back((*i)->likelihood_factor_duplicate());
    else
      this->numerator_likelihood_factors.push_back(NULL);
  }
  
  this->numerator_proposals = another.numerator_proposals;
  this->numerator_likelihood_proposals = another.numerator_likelihood_proposals;
  
  /*
  this->numerator_proposals.resize(0);
  this->numerator_proposals.reserve(another.numerator_proposals.size());
  for (auto i=another.numerator_proposals.begin();
       i!=numerator_proposals.end();
       ++i)
  {
    if (*i!=NULL)
      this->numerator_proposals.push_back((*i)->independent_proposal_kernel_duplicate());
    else
      this->numerator_proposals.push_back(NULL);
  }
  */
  //if (this->output!=NULL)
  //  this->output = another.output->duplicate();
}

// double ExactLikelihoodEstimator::estimate_log_likelihood(const List &inputs,
//                                                          const List &auxiliary_variables) const
// {
//   return this->func(inputs,this->observed_data);
// }

LikelihoodEstimatorOutput* ExactLikelihoodEstimator::initialise()
{
  return new ExactLikelihoodEstimatorOutput(this);
}

LikelihoodEstimatorOutput* ExactLikelihoodEstimator::initialise(const Parameters &parameters)
{
  return new ExactLikelihoodEstimatorOutput(this);
}

void ExactLikelihoodEstimator::setup()
{
  
}

void ExactLikelihoodEstimator::setup(const Parameters &parameters)
{
  
}

double ExactLikelihoodEstimator::evaluate(const Parameters &parameters) const
{
  double result = 0.0;

  for (auto i=this->numerator_distribution_factors.begin();
       i!=this->numerator_distribution_factors.end();
       ++i)
  {
    if (result!=-arma::datum::inf)
    {
      result = result + (*i)->evaluate(parameters);
    }
  }

  for (auto i=this->numerator_likelihood_factors.begin();
       i!=this->numerator_likelihood_factors.end();
       ++i)
  {
    if (result!=-arma::datum::inf)
    {
      result = result + (*i)->evaluate(parameters);
    }
  }
  
  for (auto i=this->numerator_proposals.begin();
       i!=this->numerator_proposals.end();
       ++i)
  {
    if (result!=-arma::datum::inf)
    {
      result = result + (*i)->evaluate_independent_kernel(parameters);
    }
  }
  
  for (auto i=this->numerator_likelihood_proposals.begin();
       i!=this->numerator_likelihood_proposals.end();
       ++i)
  {
    if (result!=-arma::datum::inf)
    {
      result = result + (*i)->evaluate_kernel(*this->current_data,
                                              parameters);
    }
  }

  return result;
}

double ExactLikelihoodEstimator::subsample_evaluate(const Parameters &parameters) const
{
  double result = 0.0;
  
  for (auto i=this->numerator_distribution_factors.begin();
       i!=this->numerator_distribution_factors.end();
       ++i)
  {
    if (result!=-arma::datum::inf)
    {
      result = result + (*i)->evaluate(parameters);
    }
  }
  
  for (auto i=this->numerator_likelihood_factors.begin();
       i!=this->numerator_likelihood_factors.end();
       ++i)
  {
    if (result!=-arma::datum::inf)
    {
      (*i)->set_data(this->subsampler->small_data);
      result = result + this->subsampler->ratio*(*i)->evaluate(parameters);
      (*i)->set_data(this->data);
    }
  }
  
  for (auto i=this->numerator_proposals.begin();
       i!=this->numerator_proposals.end();
       ++i)
  {
    if (result!=-arma::datum::inf)
    {
      result = result + (*i)->evaluate_independent_kernel(parameters);
    }
  }
  
  return result;
}

arma::mat ExactLikelihoodEstimator::evaluate_gradient(const std::string &variable,
                                                      const Parameters &parameters) const
{
  
  arma::mat parameter = parameters[variable];
  arma::mat result(parameter.n_rows,parameter.n_cols);
  result.fill(0.0);
  
  for (auto i=this->numerator_distribution_factors.begin();
       i!=this->numerator_distribution_factors.end();
       ++i)
  {
    result = result + (*i)->evaluate_gradient(variable,
                                              parameters);
  }
  
  for (auto i=this->numerator_likelihood_factors.begin();
       i!=this->numerator_likelihood_factors.end();
       ++i)
  {
    result = result + (*i)->evaluate_gradient(variable,
                                              parameters);
  }
  
  for (auto i=this->numerator_proposals.begin();
       i!=this->numerator_proposals.end();
       ++i)
  {
    result = result + (*i)->independent_gradient_of_log(variable,
                                                        parameters);
  }
  
  return result;
}

arma::mat ExactLikelihoodEstimator::subsample_evaluate_gradient(const std::string &variable,
                                                                const Parameters &parameters) const
{
  
  arma::mat parameter = parameters[variable];
  arma::mat result(parameter.n_rows,parameter.n_cols);
  result.fill(0);
  
  for (auto i=this->numerator_distribution_factors.begin();
       i!=this->numerator_distribution_factors.end();
       ++i)
  {
    result = result + (*i)->evaluate_gradient(variable,
                                              parameters);
  }
  
  for (auto i=this->numerator_likelihood_factors.begin();
       i!=this->numerator_likelihood_factors.end();
       ++i)
  {
    (*i)->set_data(this->subsampler->small_data);
    result = result + (*i)->evaluate_gradient(variable,
                                              parameters);
    (*i)->set_data(this->data);
  }
  
  for (auto i=this->numerator_proposals.begin();
       i!=this->numerator_proposals.end();
       ++i)
  {
    result = result + (*i)->subsample_independent_gradient_of_log(variable,
                                                                  parameters);
  }
  
  return result;
}

void ExactLikelihoodEstimator::specific_change_data(Data* new_data)
{
  for (auto i=this->numerator_likelihood_factors.begin();
       i!=this->numerator_likelihood_factors.end();
       ++i)
  {
    (*i)->set_data(new_data);
  }
}

//double ExactLikelihoodEstimator::evaluate(const Parameters &parameters)
//{
//  return this->func(parameters,*this->data);
//}

// void ExactLikelihoodEstimator::is_setup_likelihood_estimator(const std::vector<List> &all_points,
//                                                              const std::vector<List> &all_auxiliary_variables)
// {
//
// }
