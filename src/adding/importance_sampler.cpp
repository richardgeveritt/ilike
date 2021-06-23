#include "importance_sampler.h"

///Default constructor, Shape and Scale parameters have no values.
ImportanceSampler::ImportanceSampler(const ModelAndAlgorithm* model_and_algorithm_in)
  :SMC(model_and_algorithm_in)
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

SMC* ImportanceSampler::duplicate(void)const
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
