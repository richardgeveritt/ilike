#include "single_index.h"

SingleIndex::SingleIndex()
  :Index()
{
}

SingleIndex::~SingleIndex()
{
}

SingleIndex::SingleIndex(const SingleIndex &another)
  :Index(another)
{
  this->make_copy(another);
}

void SingleIndex::operator=(const SingleIndex &another)
{
  if(this == &another)
    return;

  Index::operator=(another);
  this->make_copy(another);
}

void SingleIndex::make_copy(const SingleIndex &another)
{
}
