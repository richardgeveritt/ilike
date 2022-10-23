#include "data_subsetter.h"
#include "distributions.h"
#include "data.h"
#include "data_subsampler.h"

DataSubsetter::DataSubsetter()
{
}

DataSubsetter::~DataSubsetter()
{
  
}

DataSubsetter::DataSubsetter(const DataSubsetter &another)
{
  this->make_copy(another);
}

void DataSubsetter::operator=(const DataSubsetter &another)
{
  if(this == &another)
    return;

  this->make_copy(another);
}

void DataSubsetter::make_copy(const DataSubsetter &another)
{
  this->data = another.data;
  this->subsampler = another.subsampler;
}
