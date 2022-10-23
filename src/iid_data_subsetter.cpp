#include "iid_data_subsetter.h"
#include "distributions.h"
#include "data.h"
#include "data_subsampler.h"

IIDDataSubsetter::IIDDataSubsetter()
  :DataSubsetter()
{
}

IIDDataSubsetter::~IIDDataSubsetter()
{
  
}

//Copy constructor for the IIDDataSubsetter class.
IIDDataSubsetter::IIDDataSubsetter(const IIDDataSubsetter &another)
  :DataSubsetter(another)
{
  this->make_copy(another);
}

void IIDDataSubsetter::operator=(const IIDDataSubsetter &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  DataSubsetter::operator=(another);
  this->make_copy(another);
}

DataSubsetter* IIDDataSubsetter::duplicate(void)const
{
  return( new IIDDataSubsetter(*this));
}

void IIDDataSubsetter::make_copy(const IIDDataSubsetter &another)
{
}

void IIDDataSubsetter::subset(size_t num_pieces)
{
  this->subsampler->ratio = 0.0;
  for (auto i=this->data->vector_begin();
       i!=this->data->vector_end();
       ++i)
  {
    arma::colvec probabilities(i->second.n_rows);
    probabilities.fill(1.0);
    arma::uvec indices = multiple_rdis(*this->subsampler->rng,
                                       num_pieces,
                                       probabilities);
    (*this->subsampler->small_data)[i->first] = (*this->data)[i->first].rows(indices);
    
    double current_ratio = double(i->second.n_rows)/double(num_pieces);
    if (this->subsampler->ratio!=0.0)
    {
      if (this->subsampler->ratio!=current_ratio)
        throw std::runtime_error("IIDDataSubsetter::subset - not yet implemented where we have different subsetting ratios for different variables in the data.");
      this->subsampler->ratio = current_ratio;
    }
    else
      this->subsampler->ratio = current_ratio;
  }
  // multinomial better in the case where we only store pointers to data
}
