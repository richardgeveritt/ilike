//
//  packing_instructions.cpp
//  ilike_cpp
//
//  Created by Richard Everitt on 21/02/2023.
//
#include "packing_instructions_map.h"

namespace ilike
{
PackingInstructionsMap::PackingInstructionsMap()
{
  
}

PackingInstructionsMap::~PackingInstructionsMap()
{
  
}

PackingInstructionsMap::PackingInstructionsMap(const PackingInstructionsMap &another)
{
  this->make_copy(another);
}

PackingInstructionsMap& PackingInstructionsMap::operator=(const PackingInstructionsMap &another)
{
  if(this == &another)
    return *this;
  
  this->make_copy(another);
  return *this;
}

void PackingInstructionsMap::make_copy(const PackingInstructionsMap &another)
{
  this->variable_indexed_start_and_end = another.variable_indexed_start_and_end;
}

void PackingInstructionsMap::set_info(const Parameters &parameters,
                                      const std::vector<std::string> &variables)
{
  size_t counter = 0;
  for (size_t i=0; i<variables.size(); ++i)
  {
    arma::colvec current = arma::vectorise(parameters[variables[i]]);
    this->variable_indexed_start_and_end[variables[i]] = std::make_pair(counter, counter+current.n_elem-1);
    counter = counter + current.n_elem;
  }
}

std::pair<size_t,size_t>& PackingInstructionsMap::operator[](const std::string &variable)
{
  return(this->variable_indexed_start_and_end[variable]);
}

std::pair<size_t,size_t> PackingInstructionsMap::operator[](const std::string &variable) const
{
  auto found =  this->variable_indexed_start_and_end.find(variable);
  
  if (found != this->variable_indexed_start_and_end.end())
    return(found->second);
  else
    return std::pair<size_t,size_t>();
}
}
