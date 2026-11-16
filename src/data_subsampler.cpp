#include "data_subsampler.h"
#include "data_subsetter.h"
#include "distributions.h"
#include "parameters.h"
#include "iid_data_subsetter.h"

namespace ilike
{
DataSubsampler::DataSubsampler()
{
  this->small_data.reset();
}

DataSubsampler::~DataSubsampler()
{
  for (auto i=this->subsetters.begin();
       i!=this->subsetters.end();
       ++i)
  {
    if (*i!=NULL)
      (*i).reset();
  }
  
  this->small_data.reset();
}

DataSubsampler::DataSubsampler(const DataSubsampler &another)
{
  this->make_copy(another);
}

void DataSubsampler::operator=(const DataSubsampler &another)
{
  if(this == &another)
    return;
  
  for (auto i=this->subsetters.begin();
       i!=this->subsetters.end();
       ++i)
  {
    if (*i!=NULL)
      (*i).reset();
  }
  
  this->subsetters.clear();
  
  this->small_data.reset();
  
  this->make_copy(another);
}

void DataSubsampler::make_copy(const DataSubsampler &another)
{
  this->variables = another.variables;
  
  this->subsetters.resize(0);
  this->subsetters.reserve(another.subsetters.size());
  for (auto i=another.subsetters.begin();
       i!=another.subsetters.end();
       ++i)
  {
    if (*i!=NULL)
      this->subsetters.emplace_back((*i)->duplicate());
    else
      this->subsetters.emplace_back();
  }
  
  if (another.small_data!=NULL)
    this->small_data.reset(another.small_data->duplicate());
  else
    this->small_data.reset();
  
  this->ratio = another.ratio;
}

void DataSubsampler::subsample(size_t num_pieces)
{
  this->small_data.reset();
  
  this->small_data.reset(new Data());
  
  for (auto i=this->subsetters.begin();
       i!=this->subsetters.end();
       ++i)
  {
    (*i)->subset(num_pieces);
  }
}
}
