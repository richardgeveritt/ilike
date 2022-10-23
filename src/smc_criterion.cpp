#include "smc_criterion.h"

SMCCriterion::SMCCriterion()
{
}

SMCCriterion::~SMCCriterion()
{

}

SMCCriterion::SMCCriterion(double desired_criterion_in)
{
  this->desired_criterion = desired_criterion_in;
}


SMCCriterion::SMCCriterion(const SMCCriterion &another)
{
  this->make_copy(another);
}

void SMCCriterion::operator=(const SMCCriterion &another)
{
  if(this == &another)
    return;

  this->make_copy(another);
}

void SMCCriterion::make_copy(const SMCCriterion &another)
{
  this->desired_criterion = another.desired_criterion;
}
