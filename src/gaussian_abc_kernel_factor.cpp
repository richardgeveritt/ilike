#include "gaussian_abc_kernel_factor.h"
#include "utils.h"
#include "distributions.h"

namespace ilike
{

GaussianABCKernelFactor::GaussianABCKernelFactor()
:ABCKernelFactor()
{
}

GaussianABCKernelFactor::GaussianABCKernelFactor(const std::vector<std::string> &data_variables_in,
                                                 const std::string &epsilon_variable_in,
                                                 Data* data_in)
:ABCKernelFactor(data_variables_in,
                 epsilon_variable_in,
                 data_in)
{
  double n = double(this->data_colvec.n_elem);
  this->constant = -n/2.0*log(2*M_PI);
}

GaussianABCKernelFactor::GaussianABCKernelFactor(const std::vector<std::string> &data_variables_in,
                                                 const std::string &epsilon_variable_in,
                                                 const std::string &scale_variable_in,
                                                 Data* data_in)
:ABCKernelFactor(data_variables_in,
                 epsilon_variable_in,
                 scale_variable_in,
                 data_in)
{
  double n = double(this->data_colvec.n_elem);
  this->constant = -(n/2.0)*log(2.0*M_PI);
}

GaussianABCKernelFactor::~GaussianABCKernelFactor()
{
}

GaussianABCKernelFactor::GaussianABCKernelFactor(const GaussianABCKernelFactor &another)
:ABCKernelFactor(another)
{
  this->make_copy(another);
}

void GaussianABCKernelFactor::operator=(const GaussianABCKernelFactor &another)
{
  if(this == &another)
    return;
  
  ABCKernelFactor::operator=(another);
  this->make_copy(another);
}

Factor* GaussianABCKernelFactor::duplicate() const
{
  return( new GaussianABCKernelFactor(*this));
}

LikelihoodFactor* GaussianABCKernelFactor::likelihood_factor_duplicate() const
{
  return( new GaussianABCKernelFactor(*this));
}

ABCKernelFactor* GaussianABCKernelFactor::abc_kernel_factor_duplicate() const
{
  return( new GaussianABCKernelFactor(*this));
}

void GaussianABCKernelFactor::make_copy(const GaussianABCKernelFactor &another)
{
  this->constant = another.constant;
}

double GaussianABCKernelFactor::likelihood_evaluate(const Parameters &input) const
{
  // At an infinite tolerance the kernel is taken to be 1 (log 0): the ABC
  // target is then the prior predictive, with evidence 1, so an SMC can start
  // from epsilon = Inf (the first weights are the kernel at the first finite
  // tolerance).
  if (std::isinf(input[this->epsilon_variable][0]))
    return 0.0;
  double epsilon2 = pow(input[this->epsilon_variable][0],2.0);
  if (epsilon2==0.0)
    return -arma::datum::inf;
  else
  {
    arma::colvec scale;
    if (this->scale_variable!="")
      scale = input[this->scale_variable];
    else
    {
      scale = arma::colvec(this->data_colvec.n_elem);
      scale.fill(1.0);
    }
    
    double distance = arma::sum(arma::pow((input.get_colvec(this->data_variables)-this->data_colvec)/scale,2.0));
    return this->constant - 0.5*(1.0/epsilon2)*distance - 0.5*scale.n_elem*log(epsilon2) - 0.5*arma::sum(log(scale));
  }
}

void GaussianABCKernelFactor::find_distance(const Parameters &input,
                                            double &distance,
                                            double &scale_constant) const
{
  arma::colvec scale;
  if (this->scale_variable!="")
  {
    scale = input[this->scale_variable];
    scale_constant = 0.5*arma::sum(log(scale));
  }
  else
  {
    scale = arma::colvec(this->data_colvec.n_elem);
    scale.fill(1.0);
    scale_constant = 0.0;
  }

  distance = arma::sum(arma::pow((input.get_colvec(this->data_variables)-this->data_colvec)/scale,2.0));
}

double GaussianABCKernelFactor::evaluate_kernel_given_distance(const Parameters &input,
                                                               double distance,
                                                               double scale_constant) const
{
  return this->evaluate_kernel_given_distance_at(input[this->epsilon_variable][0],
                                                 distance,
                                                 scale_constant);
}

double GaussianABCKernelFactor::evaluate_kernel_given_distance_at(double epsilon,
                                                                  double distance,
                                                                  double scale_constant) const
{
  // At an infinite tolerance the kernel is taken to be 1 (log 0): the ABC
  // target is then the prior predictive, with evidence 1, so an SMC can start
  // from epsilon = Inf (the first weights are the kernel at the first finite
  // tolerance).
  if (std::isinf(epsilon))
    return 0.0;
  double epsilon2 = pow(epsilon,2.0);
  if (epsilon2==0.0)
    return -arma::datum::inf;
  else
  {
    return this->constant - 0.5*(1.0/(epsilon2))*distance - 0.5*this->data_colvec.n_elem*log(epsilon2) - scale_constant;
  }
}

void GaussianABCKernelFactor::specific_set_data()
{
  this->data_colvec = this->data->get_colvec(this->data_variables);
  this->packing_instructions.set_info(*this->data,this->data_variables);
  double n = double(this->data_colvec.n_elem);
  this->constant = -n/2.0*log(2*M_PI);
}

arma::mat GaussianABCKernelFactor::likelihood_evaluate_gradient(const std::string &variable,
                                                                const Parameters &input) const
{
  double epsilon = input[this->epsilon_variable][0];

  if (this->scale_variable!="")
    return -(1.0/(epsilon*epsilon))*(input[variable]-this->data_colvec(arma::span(this->packing_instructions[variable].first,
                                                                                  this->packing_instructions[variable].second)));
  else
  {
    arma::colvec variable_scale = input[this->scale_variable](arma::span(this->packing_instructions[variable].first,
                                                                         this->packing_instructions[variable].second),0);
    return -(1.0/(epsilon*epsilon))*variable_scale%(input[variable]-this->data_colvec(arma::span(this->packing_instructions[variable].first,
                                                                                                 this->packing_instructions[variable].second)));
  }
}

}
