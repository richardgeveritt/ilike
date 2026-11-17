#include "ensemble_factors.h"
#include "vector_index.h"

namespace ilike
{
EnsembleFactors::EnsembleFactors()
{
}

EnsembleFactors::~EnsembleFactors()
{
}

EnsembleFactors::EnsembleFactors(const EnsembleFactors &another)
{
  this->make_copy(another);
}

void EnsembleFactors::operator=(const EnsembleFactors &another)
{
  if(this == &another)
    return;
  
  this->make_copy(another);
}

void EnsembleFactors::make_copy(const EnsembleFactors &another)
{
  this->measurement_names = another.measurement_names;
}

void EnsembleFactors::set_data(size_t index)
{
  VectorIndex vector_index(index);
  this->set_data(&vector_index);
}

}
