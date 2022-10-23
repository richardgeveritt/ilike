#include "ensemble_shifter.h"

EnsembleShifter::EnsembleShifter()
{
}

EnsembleShifter::~EnsembleShifter()
{

}

EnsembleShifter::EnsembleShifter(const EnsembleShifter &another)
{
  this->make_copy(another);
}

void EnsembleShifter::operator=(const EnsembleShifter &another)
{
  if(this == &another)
    return;

  this->make_copy(another);
}

void EnsembleShifter::make_copy(const EnsembleShifter &another)
{
}
