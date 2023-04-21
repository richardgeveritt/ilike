#include "standard_mcmc_output.h"
#include "particle.h"
#include "factor_variables.h"
#include "ensemble_factor_variables.h"
#include "transform.h"

StandardMCMCOutput::StandardMCMCOutput()
  :MoveOutput()
{
}

//StandardMCMCOutput::StandardMCMCOutput(const Parameters &parameters_in)
//{
//  this->output = std::deque<Particle>(1);
//  this->output[0] = Particle(parameters_in);
//}

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

MoveOutput* StandardMCMCOutput::duplicate() const
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

void StandardMCMCOutput::write_vector_points(const std::vector<std::string> &variables,
                                             std::ofstream &file_stream,
                                             std::shared_ptr<Transform> transform) const
{
  if(file_stream.is_open())
  {
    for (auto i = this->output.begin();
         i!=this->output.end();
         ++i)
    {
      if (transform==NULL)
        file_stream << i->get_rowvec(variables);
      else
        file_stream << transform->inverse_transform(i->parameters).get_rowvec(variables);
    }
  }
}

void StandardMCMCOutput::write_any_points(const std::vector<std::string> &variables,
                                          std::ofstream &file_stream) const
{
  //if(file_stream.is_open())
  //{
  //  for (auto i = this->output.begin();
  //       i!=this->output.end();
  //       ++i)
  //  {
  //    file_stream << i->parameters. << std::endl;
  //  }
  //}
}

void StandardMCMCOutput::write_factors(const std::string &directory_name,
                                       const std::string &index) const
{
  for (auto i = this->output.begin();
       i!=this->output.end();
       ++i)
  {
    if (i->factor_variables!=NULL)
    {
      i->factor_variables->write_to_file(directory_name,
                                         index);
    }
  }
}

void StandardMCMCOutput::write_ensemble_factors(const std::string &directory_name,
                                                const std::string &index) const
{
  for (auto i = this->output.begin();
       i!=this->output.end();
       ++i)
  {
    if (i->ensemble_factor_variables!=NULL)
    {
      i->ensemble_factor_variables->write_to_file(directory_name,
                                                  index);
    }
  }
}

void StandardMCMCOutput::close_ofstreams()
{
  for (auto i = this->output.begin();
       i!=this->output.end();
       ++i)
  {
    if (i->factor_variables!=NULL)
    {
      i->factor_variables->close_ofstreams();
    }
    
    if (i->ensemble_factor_variables!=NULL)
    {
      i->ensemble_factor_variables->close_ofstreams();
    }
  }
}
