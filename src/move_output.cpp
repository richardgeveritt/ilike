#include "move_output.h"

namespace ilike
{
MoveOutput::MoveOutput()
{
}

MoveOutput::~MoveOutput()
{
  
}

MoveOutput::MoveOutput(const MoveOutput &another)
{
  this->make_copy(another);
}

void MoveOutput::operator=(const MoveOutput &another)
{
  if(this == &another)
    return;
  
  this->make_copy(another);
}

void MoveOutput::make_copy(const MoveOutput &another)
{
  //this->vector_variables = another.vector_variables;
  //this->any_variables = another.any_variables;
}
}
