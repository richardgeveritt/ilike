#include "vector_factors.h"
#include "likelihood_estimator.h"
#include "vector_factor_variables.h"
#include "index.h"
#include "likelihood_estimator_output.h"

VectorFactors::VectorFactors()
  :Factors()
{
  this->likelihood_estimators.resize(0);
  this->likelihood_estimator_temp_data.resize(0);
}

VectorFactors::VectorFactors(const std::vector<LikelihoodEstimator*> &likelihood_estimators_in)
:Factors()
{
  this->likelihood_estimators = likelihood_estimators_in;
}

VectorFactors::~VectorFactors()
{
  for (std::vector<LikelihoodEstimator*>::iterator i=this->likelihood_estimators.begin();
       i!=this->likelihood_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  
  for (std::vector<Data*>::iterator i=this->likelihood_estimator_temp_data.begin();
       i!=this->likelihood_estimator_temp_data.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
}

//Copy constructor for the VectorFactors class.
VectorFactors::VectorFactors(const VectorFactors &another)
  :Factors(another)
{
  this->make_copy(another);
}

void VectorFactors::operator=(const VectorFactors &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  for (std::vector<LikelihoodEstimator*>::iterator i=this->likelihood_estimators.begin();
       i!=this->likelihood_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  this->likelihood_estimators.clear();
  
  for (std::vector<Data*>::iterator i=this->likelihood_estimator_temp_data.begin();
       i!=this->likelihood_estimator_temp_data.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  this->likelihood_estimator_temp_data.clear();
  
  Factors::operator=(another);
  this->make_copy(another);
}

Factors* VectorFactors::duplicate() const
{
  return( new VectorFactors(*this));
}

void VectorFactors::make_copy(const VectorFactors &another)
{
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
}

void VectorFactors::set_data(const Index* index)
{
  if (this->likelihood_estimators.size()==0)
    return;
  
  arma::uvec indices = index->get_uvec();
  Data* subsetted_data = new Data();
  
  // assumption that all llhd_estimators point to the same data
  Data* all_data = this->likelihood_estimators[0]->get_data();
  
  for (auto i=all_data->vector_begin();
       i!=all_data->vector_end();
       ++i)
  {
    (*subsetted_data)[i->first] = (*all_data)[i->first].rows(indices);
  }
  if (this->likelihood_estimator_temp_data[0]!=NULL)
    delete this->likelihood_estimator_temp_data[0];
  this->likelihood_estimator_temp_data[0] = subsetted_data;
  for (auto i=this->likelihood_estimators.begin();
       i!=this->likelihood_estimators.end();
       ++i)
  {
    (*i)->change_data(subsetted_data);
  }
}

FactorVariables* VectorFactors::simulate_factor_variables(const Parameters &simulated_parameters) const
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
  
  return new VectorFactorVariables(this,
                                   outputs);
}

/*
FactorVariables* VectorFactors::simulate_factor_variables(const Parameters &simulated_parameters,
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
  
  return new VectorFactorVariables(this,
                                   outputs);
}
*/

FactorVariables* VectorFactors::subsample_simulate_factor_variables(const Parameters &simulated_parameters) const
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
  
  return new VectorFactorVariables(this,
                                   outputs);
}

/*
FactorVariables* VectorFactors::subsample_simulate_factor_variables(const Parameters &simulated_parameters,
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
  
  return new VectorFactorVariables(this,
                                   outputs);
}
*/

void VectorFactors::setup()
{
  for (auto i=this->likelihood_estimators.begin();
       i!=this->likelihood_estimators.end();
       ++i)
  {
    (*i)->setup();
  }
}

void VectorFactors::setup(const Parameters &conditioned_on_parameters)
{
  for (auto i=this->likelihood_estimators.begin();
       i!=this->likelihood_estimators.end();
       ++i)
  {
    (*i)->setup(conditioned_on_parameters);
  }
}
