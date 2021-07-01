#include "importance_sampler.h"

ImportanceSampler::ImportanceSampler()
  :SMC()
{
}

ImportanceSampler::ImportanceSampler(const ModelAndAlgorithm &model_and_algorithm_in,
                                     const Data* data_in)
  :SMC(model_and_algorithm_in, data_in)
{
}

//Copy constructor for the ImportanceSampler class.
ImportanceSampler::ImportanceSampler(const ImportanceSampler &another)
  :SMC(another)
{
  this->make_copy(another);
}

//Destructor for the ImportanceSampler class.
ImportanceSampler::~ImportanceSampler(void)
{
}

void ImportanceSampler::operator=(const ImportanceSampler &another)
{
  if(this == &another){ //if a==a
    return;
  }

  SMC::operator=(another);
  this->make_copy(another);
}

SMC* ImportanceSampler::smc_duplicate(void)const
{
  return( new ImportanceSampler(*this));
}

LikelihoodEstimator* ImportanceSampler::duplicate(void)const
{
  return( new ImportanceSampler(*this));
}

void ImportanceSampler::make_copy(const ImportanceSampler &another)
{

}

void ImportanceSampler::smc_step()
{
}

void ImportanceSampler::weight_update()
{
}
