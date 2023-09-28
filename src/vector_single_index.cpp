#include <iterator>
#include "vector_single_index.h"

VectorSingleIndex::VectorSingleIndex()
  :SingleIndex()
{
}

VectorSingleIndex::VectorSingleIndex(const std::vector<size_t> &indices_in)
:SingleIndex()
{
  this->indices = indices_in;
}

VectorSingleIndex::VectorSingleIndex(size_t start_in,
                                     size_t end_in)
:SingleIndex()
{
  if (start_in>end_in)
    Rcpp::stop("RangeSingleIndex(size_t start_in, size_t end_in): start is greater than end.");
  this->indices.reserve(end_in-start_in);
  for (size_t i=start_in; i<end_in; ++i)
  {
    this->indices.push_back(i);
  }
}

VectorSingleIndex::VectorSingleIndex(size_t single_index_in)
:SingleIndex()
{
  this->indices.push_back(single_index_in);
}

VectorSingleIndex::~VectorSingleIndex()
{
}

VectorSingleIndex::VectorSingleIndex(const VectorSingleIndex &another)
  :SingleIndex(another)
{
  this->make_copy(another);
}

void VectorSingleIndex::operator=(const VectorSingleIndex &another)
{
  if(this == &another)
    return;

  SingleIndex::operator=(another);
  this->make_copy(another);
}

Index* VectorSingleIndex::duplicate() const
{
  return( new VectorSingleIndex(*this));
}

SingleIndex* VectorSingleIndex::single_index_duplicate() const
{
  return( new VectorSingleIndex(*this));
}

void VectorSingleIndex::make_copy(const VectorSingleIndex &another)
{
  this->indices = another.indices;
}

std::vector<size_t>::const_iterator VectorSingleIndex::begin() const
{
  return this->indices.begin();
}

std::vector<size_t>::const_iterator VectorSingleIndex::end() const
{
  return this->indices.end();
}

size_t VectorSingleIndex::size() const
{
  return this->indices.size();
}

arma::uvec VectorSingleIndex::get_uvec() const
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
void VectorSingleIndex::add_index(const size_t &number)
{
  this->indices.push_back(number);
}
*/
