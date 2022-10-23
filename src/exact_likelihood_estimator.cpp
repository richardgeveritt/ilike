#include "exact_likelihood_estimator.h"
#include "exact_likelihood_estimator_output.h"
#include "independent_proposal_kernel.h"

ExactLikelihoodEstimator::ExactLikelihoodEstimator()
  :LikelihoodEstimator()
{
}

ExactLikelihoodEstimator::ExactLikelihoodEstimator(RandomNumberGenerator* rng_in,
                                                   size_t* seed_in,
                                                   Data* data_in,
                                                   EvaluateLogLikelihoodPtr llhd_in,
                                                   bool smcfixed_flag_in)
:LikelihoodEstimator(rng_in, seed_in, data_in)
{
  this->numerator_llhds.reserve(1);
  this->numerator_llhds.push_back(llhd_in);
  this->smcfixed_flag = smcfixed_flag_in;
  //this->output = new ExactLikelihoodEstimatorOutput();
}

ExactLikelihoodEstimator::ExactLikelihoodEstimator(RandomNumberGenerator* rng_in,
                                                   size_t* seed_in,
                                                   Data* data_in,
                                                   EvaluateLogDistributionPtr dist_in,
                                                   bool smcfixed_flag_in)
:LikelihoodEstimator(rng_in, seed_in, data_in)
{
  this->numerator_distributions.reserve(1);
  this->numerator_distributions.push_back(dist_in);
  this->smcfixed_flag = smcfixed_flag_in;
  //this->output = new ExactDistributionEstimatorOutput();
}

ExactLikelihoodEstimator::ExactLikelihoodEstimator(RandomNumberGenerator* rng_in,
                                                   size_t* seed_in,
                                                   Data* data_in,
                                                   IndependentProposalKernel* dist_in,
                                                   bool smcfixed_flag_in)
:LikelihoodEstimator(rng_in, seed_in, data_in)
{
  this->numerator_distributions.reserve(1);
  this->numerator_proposals.push_back(dist_in);
  this->smcfixed_flag = smcfixed_flag_in;
}

ExactLikelihoodEstimator::ExactLikelihoodEstimator(RandomNumberGenerator* rng_in,
                                                   size_t* seed_in,
                                                   Data* data_in,
                                                   EvaluateLogDistributionPtr prior_in,
                                                   EvaluateLogLikelihoodPtr llhd_in,
                                                   bool smcfixed_flag_in)
:LikelihoodEstimator(rng_in, seed_in, data_in)
{
  this->numerator_llhds.reserve(1);
  this->numerator_llhds.push_back(llhd_in);
  
  this->numerator_distributions.reserve(1);
  this->numerator_distributions.push_back(prior_in);
  
  this->smcfixed_flag = smcfixed_flag_in;
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
  
  this->numerator_llhds.clear();
  this->numerator_distributions.clear();
  this->denominator_llhds.clear();
  this->denominator_distributions.clear();
  
  this->gradient_numerator_llhds.clear();
  this->gradient_numerator_distributions.clear();
  this->gradient_denominator_llhds.clear();
  this->gradient_denominator_distributions.clear();
  this->numerator_proposals.clear();

  LikelihoodEstimator::operator=(another);
  this->make_copy(another);
}

LikelihoodEstimator* ExactLikelihoodEstimator::duplicate(void)const
{
  return( new ExactLikelihoodEstimator(*this));
}

void ExactLikelihoodEstimator::make_copy(const ExactLikelihoodEstimator &another)
{
  this->numerator_llhds = another.numerator_llhds;
  this->numerator_distributions = another.numerator_distributions;
  this->denominator_llhds = another.denominator_llhds;
  this->denominator_distributions = another.denominator_distributions;
  
  this->gradient_numerator_llhds = another.gradient_numerator_llhds;
  this->gradient_numerator_distributions = another.gradient_numerator_distributions;
  this->gradient_denominator_llhds = another.gradient_denominator_llhds;
  this->gradient_denominator_distributions = another.gradient_denominator_distributions;
  
  this->numerator_proposals = another.numerator_proposals;
  
  this->smcfixed_flag = another.smcfixed_flag;
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

double ExactLikelihoodEstimator::evaluate(const Parameters &parameters)
{
  double result = 0.0;
  
  for (std::vector<EvaluateLogLikelihoodPtr>::iterator i=this->numerator_llhds.begin();
       i!=this->numerator_llhds.end();
       ++i)
  {
    result = result + (*i)(parameters,*this->data);
  }
  
  for (std::vector<EvaluateLogDistributionPtr>::iterator i=this->numerator_distributions.begin();
       i!=this->numerator_distributions.end();
       ++i)
  {
    result = result + (*i)(parameters);
  }
  
  for (auto i=this->numerator_proposals.begin();
       i!=this->numerator_proposals.end();
       ++i)
  {
    result = result + (*i)->evaluate_independent_kernel(parameters);
  }
  
  for (std::vector<EvaluateLogLikelihoodPtr>::iterator i=this->denominator_llhds.begin();
       i!=this->denominator_llhds.end();
       ++i)
  {
    result = result - (*i)(parameters,*this->data);
  }
  
  for (std::vector<EvaluateLogDistributionPtr>::iterator i=this->denominator_distributions.begin();
       i!=this->denominator_distributions.end();
       ++i)
  {
    result = result - (*i)(parameters);
  }
  
  return result;
}

double ExactLikelihoodEstimator::subsample_evaluate(const Parameters &parameters)
{
  double result = 0.0;
  
  for (std::vector<EvaluateLogLikelihoodPtr>::iterator i=this->numerator_llhds.begin();
       i!=this->numerator_llhds.end();
       ++i)
  {
    result = result + this->subsampler->ratio*(*i)(parameters,*this->subsampler->small_data);
  }
  
  for (std::vector<EvaluateLogDistributionPtr>::iterator i=this->numerator_distributions.begin();
       i!=this->numerator_distributions.end();
       ++i)
  {
    result = result + (*i)(parameters);
  }
  
  for (auto i=this->numerator_proposals.begin();
       i!=this->numerator_proposals.end();
       ++i)
  {
    result = result + (*i)->subsample_evaluate_independent_kernel(parameters);
  }
  
  for (std::vector<EvaluateLogLikelihoodPtr>::iterator i=this->denominator_llhds.begin();
       i!=this->denominator_llhds.end();
       ++i)
  {
    result = result - this->subsampler->ratio*(*i)(parameters,*this->subsampler->small_data);
  }
  
  for (std::vector<EvaluateLogDistributionPtr>::iterator i=this->denominator_distributions.begin();
       i!=this->denominator_distributions.end();
       ++i)
  {
    result = result - (*i)(parameters);
  }
  
  return result;
}

arma::mat ExactLikelihoodEstimator::evaluate_gradient(const std::string &variable,
                                                      const Parameters &parameters)
{
  
  arma::mat parameter = parameters[variable];
  arma::mat result(parameter.n_rows,parameter.n_cols);
  result.fill(0);
  
  for (std::vector<EvaluateGradientLogLikelihoodPtr>::iterator i=this->gradient_numerator_llhds.begin();
       i!=this->gradient_numerator_llhds.end();
       ++i)
  {
    result = result + (*i)(variable,
                           parameters,
                           *this->data);
  }
  
  for (std::vector<EvaluateGradientLogDistributionPtr>::iterator i=this->gradient_numerator_distributions.begin();
       i!=this->gradient_numerator_distributions.end();
       ++i)
  {
    result = result + (*i)(variable,
                           parameters);
  }
  
  for (auto i=this->numerator_proposals.begin();
       i!=this->numerator_proposals.end();
       ++i)
  {
    result = result + (*i)->independent_gradient_of_log(variable,
                                                        parameters);
  }
  
  for (std::vector<EvaluateGradientLogLikelihoodPtr>::iterator i=this->gradient_denominator_llhds.begin();
       i!=this->gradient_denominator_llhds.end();
       ++i)
  {
    result = result - (*i)(variable,
                           parameters,*this->data);
  }
  
  for (std::vector<EvaluateGradientLogDistributionPtr>::iterator i=this->gradient_denominator_distributions.begin();
       i!=this->gradient_denominator_distributions.end();
       ++i)
  {
    result = result - (*i)(variable,
                           parameters);
  }
  
  return result;
}

arma::mat ExactLikelihoodEstimator::subsample_evaluate_gradient(const std::string &variable,
                                                                const Parameters &parameters)
{
  
  arma::mat parameter = parameters[variable];
  arma::mat result(parameter.n_rows,parameter.n_cols);
  result.fill(0);
  
  for (std::vector<EvaluateGradientLogLikelihoodPtr>::iterator i=this->gradient_numerator_llhds.begin();
       i!=this->gradient_numerator_llhds.end();
       ++i)
  {
    result = result + (*i)(variable,
                           parameters,
                           *this->subsampler->small_data);
  }
  
  for (std::vector<EvaluateGradientLogDistributionPtr>::iterator i=this->gradient_numerator_distributions.begin();
       i!=this->gradient_numerator_distributions.end();
       ++i)
  {
    result = result + (*i)(variable,
                           parameters);
  }
  
  for (auto i=this->numerator_proposals.begin();
       i!=this->numerator_proposals.end();
       ++i)
  {
    result = result + (*i)->subsample_independent_gradient_of_log(variable,
                                                                  parameters);
  }
  
  for (std::vector<EvaluateGradientLogLikelihoodPtr>::iterator i=this->gradient_denominator_llhds.begin();
       i!=this->gradient_denominator_llhds.end();
       ++i)
  {
    result = result - (*i)(variable,
                           parameters,
                           *this->subsampler->small_data);
  }
  
  for (std::vector<EvaluateGradientLogDistributionPtr>::iterator i=this->gradient_denominator_distributions.begin();
       i!=this->gradient_denominator_distributions.end();
       ++i)
  {
    result = result - (*i)(variable,
                           parameters);
  }
  
  return result;
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
