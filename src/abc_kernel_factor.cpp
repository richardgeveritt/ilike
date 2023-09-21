#include "abc_kernel_factor.h"
#include "utils.h"

namespace ilike::exact_factor
{

ABCKernelFactor::ABCKernelFactor()
:LikelihoodFactor()
{
}

ABCKernelFactor::ABCKernelFactor(const std::vector<std::string> &data_variables_in,
                                 const std::string &epsilon_variable_in,
                                 Data* data_in)
:LikelihoodFactor(data_in)
{
  this->data_variables = data_variables_in;
  this->epsilon_variable = epsilon_variable_in;
  this->scale_variable = "";
  this->data_colvec = this->data->get_colvec(this->data_variables);
  this->packing_instructions.set_info(*this->data,this->data_variables);
}

ABCKernelFactor::ABCKernelFactor(const std::vector<std::string> &data_variables_in,
                                 const std::string &epsilon_variable_in,
                                 const std::string &scale_variable_in,
                                 Data* data_in)
:LikelihoodFactor(data_in)
{
  this->data_variables = data_variables_in;
  this->epsilon_variable = epsilon_variable_in;
  this->scale_variable = scale_variable_in;
  this->data_colvec = this->data->get_colvec(this->data_variables);
  this->packing_instructions.set_info(*this->data,this->data_variables);
}

ABCKernelFactor::~ABCKernelFactor()
{
}

ABCKernelFactor::ABCKernelFactor(const ABCKernelFactor &another)
:LikelihoodFactor(another)
{
  this->make_copy(another);
}

ABCKernelFactor& ABCKernelFactor::operator=(const ABCKernelFactor &another)
{
  if(this == &another)
    return *this;
  
  LikelihoodFactor::operator=(another);
  this->make_copy(another);
  
  return *this;
}

void ABCKernelFactor::make_copy(const ABCKernelFactor &another)
{
  this->data_variables = another.data_variables;
  this->epsilon_variable = another.epsilon_variable;
  this->scale_variable = another.scale_variable;
  this->data_colvec = another.data_colvec;
  this->packing_instructions = another.packing_instructions;
}

}
