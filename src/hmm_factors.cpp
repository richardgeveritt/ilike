#include "hmm_factors.h"
#include "proposal_kernel.h"
#include "hmm_factor_variables.h"
#include "likelihood_estimator.h"
#include "likelihood_estimator_output.h"

HMMFactors::HMMFactors()
  :Factors()
{
  this->likelihood_estimators.resize(0);
  this->likelihood_estimator_temp_data.resize(0);
}

HMMFactors::HMMFactors(ProposalKernel* transition_kernel_in,
                       const std::vector<LikelihoodEstimator*> &likelihood_estimators_in)
:Factors()
{
  this->transition_kernel = transition_kernel_in;
  this->likelihood_estimators = likelihood_estimators_in;
  this->likelihood_estimator_temp_data.resize(0);
}

HMMFactors::~HMMFactors()
{
  if (this->transition_kernel!=NULL)
    delete this->transition_kernel;
  
  for (std::vector<LikelihoodEstimator*>::iterator i=this->likelihood_estimators.begin();
       i!=this->likelihood_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  
  /*
  for (std::vector<Data*>::iterator i=this->likelihood_estimator_temp_data.begin();
       i!=this->likelihood_estimator_temp_data.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  */
}

//Copy constructor for the HMMFactors class.
HMMFactors::HMMFactors(const HMMFactors &another)
  :Factors(another)
{
  this->make_copy(another);
}

void HMMFactors::operator=(const HMMFactors &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  if (this->transition_kernel!=NULL)
    delete this->transition_kernel;
  
  for (std::vector<LikelihoodEstimator*>::iterator i=this->likelihood_estimators.begin();
       i!=this->likelihood_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  this->likelihood_estimators.clear();
  
   /*
  for (std::vector<Data*>::iterator i=this->likelihood_estimator_temp_data.begin();
       i!=this->likelihood_estimator_temp_data.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  this->likelihood_estimator_temp_data.clear();
  */
   
  Factors::operator=(another);
  this->make_copy(another);
}

Factors* HMMFactors::duplicate() const
{
  return( new HMMFactors(*this));
}

void HMMFactors::make_copy(const HMMFactors &another)
{
  if (another.transition_kernel!=NULL)
    this->transition_kernel = another.transition_kernel->proposal_kernel_duplicate();
  else
    this->transition_kernel = NULL;
  
  this->smcfixed_flag = another.smcfixed_flag;
  
  this->data_time_slices = another.data_time_slices;
  
  this->likelihood_estimators.resize(0);
  this->likelihood_estimators.reserve(another.likelihood_estimators.size());
  for (std::vector<LikelihoodEstimator*>::const_iterator i=another.likelihood_estimators.begin();
       i!=another.likelihood_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      this->likelihood_estimators.push_back((*i)->duplicate());
    else
      this->likelihood_estimators.push_back(NULL);
  }
  
  this->likelihood_estimator_temp_data = another.likelihood_estimator_temp_data;
  /*
  this->likelihood_estimator_temp_data.resize(0);
  for (std::vector<Data*>::const_iterator i=another.likelihood_estimator_temp_data.begin();
       i!=another.likelihood_estimator_temp_data.end();
       ++i)
  {
    if (*i!=NULL)
      this->likelihood_estimator_temp_data.push_back((*i)->duplicate());
    else
      this->likelihood_estimator_temp_data.push_back(NULL);
  }
  */
}

void HMMFactors::set_data(const Index* index)
{
  /*
  if (index->size()==1)
  {
    for (auto i=this->likelihood_estimators.begin();
         i!=this->likelihood_estimators.end();
         ++i)
    {
      (*i)->change_data(&this->data_time_slices[*index->begin()]);
    }
  }
  */
  //else
  //{
  if (this->likelihood_estimators.size()==0)
    return;
  
  arma::uvec indices = index->get_uvec();
  
  // assumption that all llhd_estimators point to the same data
  Data* all_data = this->likelihood_estimators[0]->get_data();
  
  if (this->likelihood_estimator_temp_data.size()==0)
  {
    this->likelihood_estimator_temp_data.push_back(std::make_shared<Data>(all_data->rows(indices)));
  }
  else
  {
    this->likelihood_estimator_temp_data[0] = std::make_shared<Data>(all_data->rows(indices));
  }
  
  for (auto i=this->likelihood_estimators.begin();
       i!=this->likelihood_estimators.end();
       ++i)
  {
    (*i)->change_data(this->likelihood_estimator_temp_data[0]);
  }
  //}
}

FactorVariables* HMMFactors::simulate_factor_variables(const Parameters &simulated_parameters) const
{
  std::vector<LikelihoodEstimatorOutput*> outputs;
  outputs.reserve(this->likelihood_estimators.size());
  
  for (std::vector<LikelihoodEstimator*>::const_iterator i = this->likelihood_estimators.begin();
       i != this->likelihood_estimators.end();
       ++i)
  {
    outputs.push_back((*i)->initialise(simulated_parameters));
    outputs.back()->simulate(simulated_parameters);
    outputs.back()->write_to_file_flag = false;
  }
  
  return new HMMFactorVariables(this,
                                outputs);
}

/*
FactorVariables* HMMFactors::simulate_factor_variables(const Parameters &simulated_parameters,
                                                       const Parameters &conditioned_on_parameters)
{
  std::vector<LikelihoodEstimatorOutput*> outputs;
  outputs.reserve(this->likelihood_estimators.size());
  
  Parameters all_parameters = simulated_parameters.merge(conditioned_on_parameters);
  
  for (std::vector<LikelihoodEstimator*>::const_iterator i = this->likelihood_estimators.begin();
       i != this->likelihood_estimators.end();
       ++i)
  {
    outputs.push_back((*i)->initialise(all_parameters));
    outputs.back()->simulate(all_parameters);
    outputs.back()->write_to_file_flag = false;
  }
  
  return new HMMFactorVariables(this,
                                outputs);
}
*/

FactorVariables* HMMFactors::subsample_simulate_factor_variables(const Parameters &simulated_parameters) const
{
  std::vector<LikelihoodEstimatorOutput*> outputs;
  outputs.reserve(this->likelihood_estimators.size());
  
  for (std::vector<LikelihoodEstimator*>::const_iterator i = this->likelihood_estimators.begin();
       i != this->likelihood_estimators.end();
       ++i)
  {
    outputs.push_back((*i)->initialise(simulated_parameters));
    outputs.back()->subsample_simulate(simulated_parameters);
    outputs.back()->write_to_file_flag = false;
  }
  
  return new HMMFactorVariables(this,
                                outputs);
}

/*
FactorVariables* HMMFactors::subsample_simulate_factor_variables(const Parameters &simulated_parameters,
                                                                 const Parameters &conditioned_on_parameters)
{
  std::vector<LikelihoodEstimatorOutput*> outputs;
  outputs.reserve(this->likelihood_estimators.size());
  
  Parameters all_parameters = simulated_parameters.merge(conditioned_on_parameters);
  
  for (std::vector<LikelihoodEstimator*>::const_iterator i = this->likelihood_estimators.begin();
       i != this->likelihood_estimators.end();
       ++i)
  {
    outputs.push_back((*i)->initialise(all_parameters));
    outputs.back()->subsample_simulate(all_parameters);
    outputs.back()->write_to_file_flag = false;
  }
  
  return new HMMFactorVariables(this,
                                outputs);
}
*/

void HMMFactors::setup()
{
  for (auto i=this->likelihood_estimators.begin();
       i!=this->likelihood_estimators.end();
       ++i)
  {
    (*i)->setup();
  }
}

void HMMFactors::setup(const Parameters &conditioned_on_parameters)
{
  for (auto i=this->likelihood_estimators.begin();
       i!=this->likelihood_estimators.end();
       ++i)
  {
    (*i)->setup(conditioned_on_parameters);
  }
}
