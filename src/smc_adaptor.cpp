#include "smc_adaptor.h"

namespace ilike
{
SMCAdaptor::SMCAdaptor()
{
}

SMCAdaptor::~SMCAdaptor()
{
}

SMCAdaptor::SMCAdaptor(const SMCAdaptor &another)
{
  this->make_copy(another);
}

void SMCAdaptor::operator=(const SMCAdaptor &another)
{
  if(this == &another)
    return;
  
  this->make_copy(another);
}

void SMCAdaptor::make_copy(const SMCAdaptor &another)
{
}
}
