//
//  packing_instructions.cpp
//  ilike_cpp
//
//  Created by Richard Everitt on 21/02/2023.
//
#include "packing_instructions.h"

namespace ilike
{
PackingInstructions::PackingInstructions()
{
  
}

PackingInstructions::~PackingInstructions()
{
  
}

PackingInstructions::PackingInstructions(const PackingInstructions &another)
{
  this->make_copy(another);
}

PackingInstructions& PackingInstructions::operator=(const PackingInstructions &another)
{
  if(this == &another)
    return *this;
  
  this->make_copy(another);
  return *this;
}

void PackingInstructions::make_copy(const PackingInstructions &another)
{
  this->states_names = another.states_names;
  this->states_start_and_end = another.states_start_and_end;
}

void PackingInstructions::set_info(const Parameters &parameters,
                                   const std::vector<std::string> &variables)
{
  this->states_names = variables;
  this->states_start_and_end.clear();
  this->states_start_and_end.reserve(variables.size());
  size_t counter = 0;
  for (size_t i=0; i<variables.size(); ++i)
  {
    arma::colvec current = arma::vectorise(parameters[variables[i]]);
    this->states_start_and_end.push_back(std::make_pair(counter, counter+current.n_elem-1));
    counter = counter + current.n_elem;
  }
}

}
