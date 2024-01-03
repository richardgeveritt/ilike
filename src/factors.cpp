#include "factors.h"
#include "vector_index.h"

Factors::Factors()
{
}

Factors::~Factors()
{
}

Factors::Factors(const Factors &another)
{
  this->make_copy(another);
}

void Factors::operator=(const Factors &another)
{
  if(this == &another)
    return;

  this->make_copy(another);
}

void Factors::make_copy(const Factors &another)
{
}

void Factors::set_data(size_t index)
{
  VectorIndex vector_index(index);
  this->set_data(&vector_index);
}

void Factors::change_data(Data* new_data)
{
  this->specific_change_data(new_data);
}
