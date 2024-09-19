#include "scale.h"

namespace ilike
{
Scale::Scale()
{
}

Scale::~Scale()
{
}

Scale::Scale(double constant_in)
{
  this->constant = constant_in;
  this->divide_by_dimension = false;
  this->dimension = 1.0;
}

Scale::Scale(double constant_in,
             size_t dimension_in)
{
  this->constant = constant_in;
  this->divide_by_dimension = true;
  this->dimension = double(dimension_in);
}

Scale::Scale(const Scale &another)
{
  this->make_copy(another);
}

Scale& Scale::operator=(const Scale &another)
{
  if(this == &another)
    return *this;
  
  this->make_copy(another);
  return *this;
}

void Scale::make_copy(const Scale &another)
{
  this->constant = another.constant;
  this->divide_by_dimension = another.divide_by_dimension;
  this->dimension = another.dimension;
}

Scale::Scale(Scale &&another)
{
  this->make_copy(std::move(another));
}

Scale& Scale::operator=(Scale &&another)
{
  if(this == &another)
    return *this;
  
  this->make_copy(std::move(another));
  return *this;
}

void Scale::make_copy(Scale &&another)
{
  this->constant = std::move(another.constant);
  this->divide_by_dimension = std::move(another.divide_by_dimension);
  this->dimension = std::move(another.dimension);
  
  another.constant = 0.0;
  another.divide_by_dimension = false;
  another.dimension = 0.0;
}

double Scale::operator()() const
{
  if (this->divide_by_dimension)
    return this->constant/this->dimension;
  else
    return this->constant;
}

double Scale::get_constant() const
{
  return this->constant;
}

double& Scale::get_constant()
{
  return this->constant;
}

/*
 double Scale::operator()(size_t dimension)
 {
 if (this->divide_by_dimension)
 return this->constant/double(dimension);
 else
 return this->constant;
 }
 */
}
