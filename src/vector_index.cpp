#include <iterator>
#include "vector_index.h"

VectorIndex::VectorIndex()
  :Index()
{
}

VectorIndex::VectorIndex(const std::vector<size_t> &indices_in)
:Index()
{
  this->indices = indices_in;
}

VectorIndex::VectorIndex(size_t start_in,
                                     size_t end_in)
:Index()
{
  if (start_in>end_in)
    Rcpp::stop("RangeIndex(size_t start_in, size_t end_in): start is greater than end.");
  this->indices.reserve(end_in-start_in);
  for (size_t i=start_in; i<end_in; ++i)
  {
    this->indices.push_back(i);
  }
}

VectorIndex::VectorIndex(size_t single_index_in)
:Index()
{
  this->indices.push_back(single_index_in);
}

VectorIndex::~VectorIndex()
{
}

VectorIndex::VectorIndex(const VectorIndex &another)
  :Index(another)
{
  this->make_copy(another);
}

void VectorIndex::operator=(const VectorIndex &another)
{
  if(this == &another)
    return;

  Index::operator=(another);
  this->make_copy(another);
}

Index* VectorIndex::duplicate() const
{
  return( new VectorIndex(*this));
}

Index* VectorIndex::single_index_duplicate() const
{
  return( new VectorIndex(*this));
}

std::vector<size_t> VectorIndex::get_indices() const
{
  return this->indices;
}

bool VectorIndex::get_transition_model() const
{
  return false;
}

void VectorIndex::make_copy(const VectorIndex &another)
{
  this->indices = another.indices;
}

std::vector<size_t>::const_iterator VectorIndex::begin() const
{
  return this->indices.begin();
}

std::vector<size_t>::const_iterator VectorIndex::end() const
{
  return this->indices.end();
}

size_t VectorIndex::size() const
{
  return this->indices.size();
}

arma::uvec VectorIndex::get_uvec() const
{
  arma::uvec result(this->indices.size());
  for (size_t i=0;
       i<this->indices.size();
       ++i)
  {
    result[i] = this->indices[i];
  }
  return result;
}

/*
void VectorIndex::add_index(const size_t &number)
{
  this->indices.push_back(number);
}
*/
