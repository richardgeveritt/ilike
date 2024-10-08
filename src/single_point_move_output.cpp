#include "single_point_move_output.h"
#include "particle.h"
#include "factor_variables.h"
#include "ensemble_factor_variables.h"
#include "ensemble_factors.h"
#include "transform.h"

namespace ilike
{
SinglePointMoveOutput::SinglePointMoveOutput()
:MoveOutput()
{
}

SinglePointMoveOutput::SinglePointMoveOutput(const Parameters &parameters_in,
                                             const Factors* factors_in,
                                             const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
                                             const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in)
: output(parameters_in,factors_in,proposals_to_transform_for_in,proposals_to_find_gradient_for_in)
{
}

SinglePointMoveOutput::SinglePointMoveOutput(const Parameters &parameters_in,
                                             const EnsembleFactors* factors_in)
: output(parameters_in,factors_in)
{
  
}

SinglePointMoveOutput::SinglePointMoveOutput(const Particle &particle_in)
{
  this->output = particle_in;
}

SinglePointMoveOutput::SinglePointMoveOutput(Parameters &&parameters_in,
                                             const Factors* factors_in,
                                             const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
                                             const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in)
: output(std::move(parameters_in),factors_in,proposals_to_transform_for_in,proposals_to_find_gradient_for_in)
{
}

SinglePointMoveOutput::SinglePointMoveOutput(Parameters &&parameters_in,
                                             const EnsembleFactors* factors_in)
: output(std::move(parameters_in),factors_in)
{
  
}

SinglePointMoveOutput::SinglePointMoveOutput(Particle &&particle_in)
{
  this->output = std::move(particle_in);
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

MoveOutput* SinglePointMoveOutput::duplicate() const
{
  return( new SinglePointMoveOutput(*this));
}

void SinglePointMoveOutput::make_copy(const SinglePointMoveOutput &another)
{
  this->output = another.output;
  this->algorithm_parameters = another.algorithm_parameters;
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
  local_output.push_back(this->output.parameters);
  return local_output;
}

void SinglePointMoveOutput::write_vector_points(const std::vector<std::string> &variables,
                                                std::ofstream &file_stream,
                                                std::shared_ptr<Transform> transform) const
{
  if (file_stream.is_open())
  {
    if (transform==NULL)
    {
      //for (size_t i=0; i<variables.size(); ++i)
      //  file_stream << output.get_rowvec(variables[i]) << ";";
      file_stream.precision(std::numeric_limits<double>::max_digits10);
      output.get_rowvec(variables).raw_print(file_stream);
      
      //file_stream << std::fixed << std::setprecision(12) << output.get_rowvec(variables);
    }
    else
    {
      //for (size_t i=0; i<variables.size(); ++i)
      //  file_stream << inverse_transform(output.parameters).get_rowvec(variables[i]) << ";";
      
      file_stream.precision(std::numeric_limits<double>::max_digits10);
      transform->inverse_transform(output.parameters).get_rowvec(variables).raw_print(file_stream);
      
      //file_stream << std::fixed << std::setprecision(12) << transform->inverse_transform(output.parameters).get_rowvec(variables);
    }
  }
}

void SinglePointMoveOutput::write_any_points(const std::vector<std::string> &variables,
                                             std::ofstream &file_stream) const
{
  
}

void SinglePointMoveOutput::write_factors(const std::string &directory_name,
                                          const std::string &index) const
{
  if (output.factor_variables!=NULL)
  {
    output.factor_variables->write_to_file(directory_name,
                                           index);
  }
}

void SinglePointMoveOutput::write_ensemble_factors(const std::string &directory_name,
                                                   const std::string &index) const
{
  if (output.ensemble_factor_variables!=NULL)
  {
    output.ensemble_factor_variables->write_to_file(directory_name,
                                                    index);
  }
}

size_t SinglePointMoveOutput::length() const
{
  return 1;
}

void SinglePointMoveOutput::close_ofstreams()
{
  if (output.factor_variables!=NULL)
  {
    output.factor_variables->close_ofstreams();
  }
  
  if (output.ensemble_factor_variables!=NULL)
  {
    output.ensemble_factor_variables->close_ofstreams();
  }
}

Parameters SinglePointMoveOutput::get_current_algorithm_parameters() const
{
  return this->algorithm_parameters;
}
}

