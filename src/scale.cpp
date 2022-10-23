#include "scale.h"

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
  
void Scale::operator=(const Scale &another)
{
  if(this == &another)
    return;

  this->make_copy(another);
}

void Scale::make_copy(const Scale &another)
{
  this->constant = another.constant;
  this->divide_by_dimension = another.divide_by_dimension;
  this->dimension = another.dimension;
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
