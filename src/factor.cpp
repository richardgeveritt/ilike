#include "factor.h"

namespace ilike
{
Factor::Factor()
{
}

Factor::~Factor()
{
}

Factor::Factor(const Factor &another)
{
  this->make_copy(another);
}

void Factor::operator=(const Factor &another)
{
  if(this == &another)
    return;
  
  this->make_copy(another);
}

void Factor::make_copy(const Factor &another)
{
  
}
}
