#include "index.h"

Index::Index()
{
}

Index::~Index()
{
}

Index::Index(const Index &another)
{
  this->make_copy(another);
}

void Index::operator=(const Index &another)
{
  if(this == &another)
    return;

  this->make_copy(another);
}

void Index::make_copy(const Index &another)
{
  
}
