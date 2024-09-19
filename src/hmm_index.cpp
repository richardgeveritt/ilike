#include <iterator>
#include "hmm_index.h"

namespace ilike
{
HMMIndex::HMMIndex()
:Index()
{
}

HMMIndex::HMMIndex(const std::vector<size_t> &likelihood_indices_time_zero_in,
                   //bool evaluate_prior_in,
                   const std::vector<size_t> &likelihood_indices_time_not_zero_in,
                   bool evaluate_transition_model_in,
                   size_t first_time_index_in)
:Index()
{
  this->likelihood_indices_time_zero = likelihood_indices_time_zero_in;
  //this->evaluate_prior = evaluate_prior;
  this->likelihood_indices_time_not_zero = likelihood_indices_time_not_zero_in;
  
  this->evaluate_transition_model = evaluate_transition_model_in;
  this->first_time_index = first_time_index_in;
}

HMMIndex::~HMMIndex()
{
}

HMMIndex::HMMIndex(const HMMIndex &another)
:Index(another)
{
  this->make_copy(another);
}

void HMMIndex::operator=(const HMMIndex &another)
{
  if(this == &another)
    return;
  
  Index::operator=(another);
  this->make_copy(another);
}

Index* HMMIndex::duplicate() const
{
  return( new HMMIndex(*this));
}

HMMIndex* HMMIndex::hmm_index_duplicate() const
{
  return( new HMMIndex(*this));
}

void HMMIndex::make_copy(const HMMIndex &another)
{
  this->time_index = another.time_index;
  this->first_time_index = another.first_time_index;
  this->likelihood_indices_time_zero = another.likelihood_indices_time_zero;
  //this->evaluate_prior = another.evaluate_prior;
  this->likelihood_indices_time_not_zero = another.likelihood_indices_time_not_zero;
  this->evaluate_transition_model = another.evaluate_transition_model;
}

void HMMIndex::set_time_index(size_t time_index_in)
{
  this->time_index = time_index_in;
}

std::vector<size_t>::const_iterator HMMIndex::begin() const
{
  if (this->time_index==this->first_time_index)
    return this->likelihood_indices_time_zero.begin();
  else
    return this->likelihood_indices_time_not_zero.begin();
}

std::vector<size_t>::const_iterator HMMIndex::end() const
{
  if (this->time_index==this->first_time_index)
    return this->likelihood_indices_time_zero.end();
  else
    return this->likelihood_indices_time_not_zero.end();
}

size_t HMMIndex::size() const
{
  if (this->time_index==this->first_time_index)
    return this->likelihood_indices_time_zero.size();
  else
    return this->likelihood_indices_time_not_zero.size();
}

arma::uvec HMMIndex::get_uvec() const
{
  
  if (this->time_index==this->first_time_index)
  {
    arma::uvec result(this->likelihood_indices_time_zero.size());
    
    for (size_t i=0;
         i<this->likelihood_indices_time_zero.size();
         ++i)
    {
      result[i] = this->likelihood_indices_time_zero[i];
    }
    
    return result;
    
  }
  else
  {
    arma::uvec result(this->likelihood_indices_time_not_zero.size());
    
    for (size_t i=0;
         i<this->likelihood_indices_time_not_zero.size();
         ++i)
    {
      result[i] = this->likelihood_indices_time_not_zero[i];
    }
    return result;
  }
  
}

bool HMMIndex::get_transition_model() const
{
  if (this->time_index==this->first_time_index)
    return false;
  else
    return this->evaluate_transition_model;
}

/*
 std::vector<size_t> HMMIndex::get_prior_indices() const
 {
 return this->prior_indices;
 }
 
 std::vector<size_t> HMMIndex::get_transition_indices() const
 {
 return this->transition_indices;
 }
 
 std::vector<size_t> HMMIndex::get_likelihood_indices() const
 {
 return this->likelihood_indices;
 }
 */

/*
 void HMMIndex::add_index(const size_t &number)
 {
 this->indices.push_back(number);
 }
 */
}
