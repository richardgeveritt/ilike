#include "data_subsampler.h"
#include "data_subsetter.h"
#include "distributions.h"
#include "parameters.h"
#include "iid_data_subsetter.h"

DataSubsampler::DataSubsampler()
{
  this->small_data = NULL;
}

DataSubsampler::~DataSubsampler()
{
  for (std::vector<DataSubsetter*>::iterator i=this->subsetters.begin();
       i!=this->subsetters.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  
  if (this->small_data!=NULL)
    delete this->small_data;
}

DataSubsampler::DataSubsampler(const DataSubsampler &another)
{
  this->make_copy(another);
}

void DataSubsampler::operator=(const DataSubsampler &another)
{
  if(this == &another)
    return;
  
  for (std::vector<DataSubsetter*>::iterator i=this->subsetters.begin();
       i!=this->subsetters.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  
  this->subsetters.clear();
  
  if (this->small_data!=NULL)
    delete this->small_data;

  this->make_copy(another);
}

void DataSubsampler::make_copy(const DataSubsampler &another)
{
  this->variables = another.variables;
  
  this->subsetters.resize(0);
  this->subsetters.reserve(another.subsetters.size());
  for (std::vector<DataSubsetter*>::const_iterator i=another.subsetters.begin();
       i!=another.subsetters.end();
       ++i)
  {
    if (*i!=NULL)
      this->subsetters.push_back((*i)->duplicate());
    else
      this->subsetters.push_back(NULL);
  }
  
  if (another.small_data!=NULL)
    this->small_data = another.small_data->duplicate();
  else
    this->small_data = NULL;
  
  this->ratio = another.ratio;
}

void DataSubsampler::subsample(size_t num_pieces)
{
  if (this->small_data!=NULL)
    delete small_data;
  
  this->small_data = new Data();
  
  for (auto i=this->subsetters.begin();
       i!=this->subsetters.end();
       ++i)
  {
    (*i)->subset(num_pieces);
  }
}
