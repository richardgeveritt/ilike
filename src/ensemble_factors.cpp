#include "ensemble_factors.h"
#include "vector_single_index.h"

EnsembleFactors::EnsembleFactors()
{
  this->temperature = 1.0;
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
  this->temperature = another.temperature;
}

void EnsembleFactors::set_data(size_t index)
{
  VectorSingleIndex vector_single_index(index);
  this->set_data(&vector_single_index);
}

void EnsembleFactors::set_temperature(double temperature_in)
{
  this->temperature = temperature_in;
}

double EnsembleFactors::get_temperature() const
{
  return this->temperature;
}
