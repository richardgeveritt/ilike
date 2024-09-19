#include "kernel.h"

namespace ilike
{
Kernel::Kernel()
{
}

Kernel::~Kernel()
{
}

Kernel::Kernel(const Kernel &another)
{
  this->make_copy(another);
}

void Kernel::operator=(const Kernel &another)
{
  if(this == &another)
    return;
  
  this->make_copy(another);
}

void Kernel::make_copy(const Kernel &another)
{
}
}
