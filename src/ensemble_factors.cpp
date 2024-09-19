#include "ensemble_factors.h"
#include "vector_index.h"

namespace ilike
{
EnsembleFactors::EnsembleFactors()
{
  //this->temperature = 0.0;
  //this->previous_temperature = 0.0;
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
  //this->temperature = another.temperature;
  //this->previous_temperature = another.previous_temperature;
}

void EnsembleFactors::set_data(size_t index)
{
  VectorIndex vector_index(index);
  this->set_data(&vector_index);
}

/*
 void EnsembleFactors::set_temperature(double temperature_in)
 {
 this->previous_temperature = this->temperature;
 this->temperature = temperature_in;
 }
 
 double EnsembleFactors::get_temperature() const
 {
 return this->temperature;
 }
 
 double EnsembleFactors::get_inverse_incremental_temperature() const
 {
 return 1.0/(this->temperature - this->previous_temperature);
 }
 */
}
