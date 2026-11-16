#include "proposal_store.h"
#include "gradient_estimator_output.h"

namespace ilike
{
ProposalStore::ProposalStore()
{
  this->gradient_estimator_output.reset();
}

ProposalStore::~ProposalStore()
{
  this->gradient_estimator_output.reset();
}

ProposalStore::ProposalStore(const Parameters &transformed_parameters_in)
{
  this->transformed_parameters = transformed_parameters_in;
  this->gradient_estimator_output.reset();
}

ProposalStore::ProposalStore(GradientEstimatorOutput* gradient_estimator_output_in)
{
  this->gradient_estimator_output.reset(gradient_estimator_output_in);
}

ProposalStore::ProposalStore(const Parameters &transformed_parameters_in,
                             GradientEstimatorOutput* gradient_estimator_output_in)
{
  this->transformed_parameters = transformed_parameters_in;
  this->gradient_estimator_output.reset(gradient_estimator_output_in);
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
  
  this->gradient_estimator_output.reset();
  
  this->make_copy(another);
}

ProposalStore::ProposalStore(ProposalStore &&another)
{
  this->transformed_parameters = std::move(another.transformed_parameters);
  this->gradient_estimator_output = std::move(another.gradient_estimator_output);
}

void ProposalStore::operator=(ProposalStore &&another)
{
  if(this == &another)
    return;
  this->gradient_estimator_output.reset();
  this->transformed_parameters = std::move(another.transformed_parameters);
  this->gradient_estimator_output = std::move(another.gradient_estimator_output);
}

void ProposalStore::make_copy(const ProposalStore &another)
{
  this->transformed_parameters = another.transformed_parameters;
  
  if (another.gradient_estimator_output!=NULL)
  {
    this->gradient_estimator_output.reset(another.gradient_estimator_output->duplicate());
  }
  else
  {
    this->gradient_estimator_output.reset();
  }
}

void ProposalStore::set_transformed_parameters(const Parameters &transformed_parameters_in)
{
  this->transformed_parameters = transformed_parameters_in;
}

void ProposalStore::set_gradient_estimator_output(GradientEstimatorOutput* gradient_estimator_output_in)
{
  this->gradient_estimator_output.reset(gradient_estimator_output_in);
}

GradientEstimatorOutput* ProposalStore::get_gradient_estimator_output() const
{
  return this->gradient_estimator_output.get();
}

const Parameters& ProposalStore::get_transformed_parameters() const
{
  return this->transformed_parameters;
}
}
