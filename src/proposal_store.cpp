#include "proposal_store.h"
#include "gradient_estimator_output.h"

ProposalStore::ProposalStore()
{
  this->gradient_estimator_output = NULL;
}

ProposalStore::~ProposalStore()
{
  if (this->gradient_estimator_output!=NULL)
    delete this->gradient_estimator_output;
}

ProposalStore::ProposalStore(const Parameters &transformed_parameters_in)
{
  this->transformed_parameters = transformed_parameters_in;
  this->gradient_estimator_output = NULL;
}

ProposalStore::ProposalStore(GradientEstimatorOutput* gradient_estimator_output_in)
{
  this->gradient_estimator_output = gradient_estimator_output_in;
}

ProposalStore::ProposalStore(const Parameters &transformed_parameters_in,
                             GradientEstimatorOutput* gradient_estimator_output_in)
{
  this->transformed_parameters = transformed_parameters_in;
  this->gradient_estimator_output = gradient_estimator_output_in;
}

ProposalStore::ProposalStore(const ProposalStore &another)
{
  this->make_copy(another);
}

void ProposalStore::operator=(const ProposalStore &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  if (this->gradient_estimator_output!=NULL)
    delete this->gradient_estimator_output;
  
  this->make_copy(another);
}

void ProposalStore::make_copy(const ProposalStore &another)
{
  this->transformed_parameters = another.transformed_parameters;
  
  if (another.gradient_estimator_output!=NULL)
  {
    this->gradient_estimator_output = another.gradient_estimator_output->duplicate();
  }
  else
  {
    this->gradient_estimator_output = NULL;
  }
}

void ProposalStore::set_transformed_parameters(const Parameters &transformed_parameters_in)
{
  this->transformed_parameters = transformed_parameters_in;
}

void ProposalStore::set_gradient_estimator_output(GradientEstimatorOutput* gradient_estimator_output_in)
{
  this->gradient_estimator_output = gradient_estimator_output_in;
}

GradientEstimatorOutput* ProposalStore::get_gradient_estimator_output() const
{
  return this->gradient_estimator_output;
}

Parameters ProposalStore::get_transformed_parameters() const
{
  return this->transformed_parameters;
}
