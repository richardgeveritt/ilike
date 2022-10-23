#include "single_point_move_output.h"
#include "particle.h"

SinglePointMoveOutput::SinglePointMoveOutput()
  :MoveOutput()
{
}

SinglePointMoveOutput::SinglePointMoveOutput(const Parameters &parameters_in)
{
  this->output = Particle(parameters_in);
}

SinglePointMoveOutput::SinglePointMoveOutput(const Particle &particle_in)
{
  this->output = particle_in;
}

SinglePointMoveOutput::~SinglePointMoveOutput()
{
  
}

//Copy constructor for the SinglePointMoveOutput class.
SinglePointMoveOutput::SinglePointMoveOutput(const SinglePointMoveOutput &another)
  :MoveOutput(another)
{
  this->make_copy(another);
}

void SinglePointMoveOutput::operator=(const SinglePointMoveOutput &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  MoveOutput::operator=(another);
  this->make_copy(another);
}

MoveOutput* SinglePointMoveOutput::duplicate(void)const
{
  return( new SinglePointMoveOutput(*this));
}

void SinglePointMoveOutput::make_copy(const SinglePointMoveOutput &another)
{
  this->output = another.output;
}

Particle& SinglePointMoveOutput::back()
{
  return this->output;
}

Particle SinglePointMoveOutput::back() const
{
  return this->output;
}

std::vector<Parameters> SinglePointMoveOutput::get_vector_of_parameters() const
{
  std::vector<Parameters> local_output;
  local_output.reserve(1);
  local_output.push_back(*this->output.move_parameters);
  return local_output;
}
