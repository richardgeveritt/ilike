#include "standard_mcmc_output.h"
#include "particle.h"

StandardMCMCOutput::StandardMCMCOutput()
  :MoveOutput()
{
}

StandardMCMCOutput::StandardMCMCOutput(const Parameters &parameters_in)
{
  this->output = std::deque<Particle>(1);
  this->output[0] = Particle(parameters_in);
}

StandardMCMCOutput::~StandardMCMCOutput()
{
  
}

//Copy constructor for the StandardMCMCOutput class.
StandardMCMCOutput::StandardMCMCOutput(const StandardMCMCOutput &another)
  :MoveOutput(another)
{
  this->make_copy(another);
}

void StandardMCMCOutput::operator=(const StandardMCMCOutput &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  MoveOutput::operator=(another);
  this->make_copy(another);
}

MoveOutput* StandardMCMCOutput::duplicate(void)const
{
  return( new StandardMCMCOutput(*this));
}

void StandardMCMCOutput::make_copy(const StandardMCMCOutput &another)
{
  this->output = another.output;
}

void StandardMCMCOutput::push_back(const Particle &particle_in)
{
  this->output.push_back(particle_in);
}

Particle& StandardMCMCOutput::back()
{
  return this->output.back();
}

Particle StandardMCMCOutput::back() const
{
  return this->output.back();
}

std::vector<Parameters> StandardMCMCOutput::get_vector_of_parameters() const
{
  std::vector<Parameters> local_output;
  local_output.reserve(this->output.size());
  for (auto i=this->output.begin(); i!=this->output.end(); ++i)
  {
    local_output.push_back(*i->move_parameters);
  }
  return local_output;
}
