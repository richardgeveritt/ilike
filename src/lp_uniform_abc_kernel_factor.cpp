#include "lp_uniform_abc_kernel_factor.h"
#include "utils.h"

namespace ilike
{

// -log of the volume of the unit ball {x in R^n : ||x||_p <= 1}:
// (2 Gamma(1/p + 1))^n / Gamma(n/p + 1); for p = Inf the cube [-1, 1]^n
static double lp_unit_ball_log_inverse_volume(double n, double p)
{
  if (std::isinf(p))
    return -n*log(2.0);
  return lgamma(n/p + 1.0) - n*(log(2.0) + lgamma(1/p + 1));
}


LpUniformABCKernelFactor::LpUniformABCKernelFactor()
:ABCKernelFactor()
{
}

LpUniformABCKernelFactor::LpUniformABCKernelFactor(double p_in,
                                                   const std::vector<std::string> &data_variables_in,
                                                   const std::string &epsilon_variable_in,
                                                   Data* data_in)
:ABCKernelFactor(data_variables_in,
                 epsilon_variable_in,
                 data_in)
{
  this->p = p_in;
  double n = double(this->data_colvec.n_elem);
  this->constant = lp_unit_ball_log_inverse_volume(n, this->p);
}

LpUniformABCKernelFactor::LpUniformABCKernelFactor(double p_in,
                                                   const std::vector<std::string> &data_variables_in,
                                                   const std::string &epsilon_variable_in,
                                                   const std::string &scale_variable_in,
                                                   Data* data_in)
:ABCKernelFactor(data_variables_in,
                 epsilon_variable_in,
                 scale_variable_in,
                 data_in)
{
  this->p = p_in;
  double n = double(this->data_colvec.n_elem);
  this->constant = lp_unit_ball_log_inverse_volume(n, this->p);
}

LpUniformABCKernelFactor::~LpUniformABCKernelFactor()
{
}

LpUniformABCKernelFactor::LpUniformABCKernelFactor(const LpUniformABCKernelFactor &another)
:ABCKernelFactor(another)
{
  this->make_copy(another);
}

void LpUniformABCKernelFactor::operator=(const LpUniformABCKernelFactor &another)
{
  if(this == &another)
    return;
  
  ABCKernelFactor::operator=(another);
  this->make_copy(another);
}

Factor* LpUniformABCKernelFactor::duplicate() const
{
  return( new LpUniformABCKernelFactor(*this));
}

LikelihoodFactor* LpUniformABCKernelFactor::likelihood_factor_duplicate() const
{
  return( new LpUniformABCKernelFactor(*this));
}

ABCKernelFactor* LpUniformABCKernelFactor::abc_kernel_factor_duplicate() const
{
  return( new LpUniformABCKernelFactor(*this));
}

void LpUniformABCKernelFactor::make_copy(const LpUniformABCKernelFactor &another)
{
  this->data_variables = another.data_variables;
  this->epsilon_variable = another.epsilon_variable;
  this->data_colvec = another.data_colvec;
  this->constant = another.constant;
  this->p = another.p;
}

double LpUniformABCKernelFactor::likelihood_evaluate(const Parameters &input) const
{
  // At an infinite tolerance the kernel is taken to be 1 (log 0): the ABC
  // target is then the prior predictive, with evidence 1, so an SMC can start
  // from epsilon = Inf (the first weights are the kernel at the first finite
  // tolerance).
  if (std::isinf(input[this->epsilon_variable][0]))
    return 0.0;
  arma::colvec scale;
  if (this->scale_variable!="")
    scale = input[this->scale_variable];
  else
  {
    scale = arma::colvec(this->data_colvec.n_elem);
    scale.fill(1.0);
  }
  double distance = scaled_Lp_distance(input.get_colvec(this->data_variables),
                                       this->data_colvec,
                                       scale,
                                       this->p);
  
  double epsilon = input[this->epsilon_variable][0];
  if (distance>epsilon)
    return -arma::datum::inf;
  else
  {
    return this->constant - scale.n_elem*log(epsilon) - arma::sum(log(scale));
  }
}

void LpUniformABCKernelFactor::find_distance(const Parameters &input,
                                             double &distance,
                                             double &scale_constant) const
{
  arma::colvec scale;
  if (this->scale_variable!="")
  {
    scale = input[this->scale_variable];
    scale_constant = arma::sum(log(scale));
  }
  else
  {
    scale = arma::colvec(this->data_colvec.n_elem);
    scale.fill(1.0);
    scale_constant = 0.0;
  }
  distance = scaled_Lp_distance(input.get_colvec(this->data_variables),
                                this->data_colvec,
                                scale,
                                this->p);
}

double LpUniformABCKernelFactor::evaluate_kernel_given_distance(const Parameters &input,
                                                                double distance,
                                                                double scale_constant) const
{
  // At an infinite tolerance the kernel is taken to be 1 (log 0): the ABC
  // target is then the prior predictive, with evidence 1, so an SMC can start
  // from epsilon = Inf (the first weights are the kernel at the first finite
  // tolerance).
  return this->evaluate_kernel_given_distance_at(input[this->epsilon_variable][0],
                                                 distance,
                                                 scale_constant);
}

double LpUniformABCKernelFactor::evaluate_kernel_given_distance_at(double epsilon,
                                                                   double distance,
                                                                   double scale_constant) const
{
  if (std::isinf(epsilon))
    return 0.0;
  if (distance>epsilon)
    return -arma::datum::inf;
  else
  {
    return this->constant - this->data_colvec.n_elem*log(epsilon) - scale_constant;
  }
}

void LpUniformABCKernelFactor::specific_set_data()
{
  this->data_colvec = this->data->get_colvec(this->data_variables);
  this->packing_instructions.set_info(*this->data,this->data_variables);
  double n = double(this->data_colvec.n_elem);
  // (the unit Lp ball's volume for this p: it used to be the Euclidean
  // ball's here, right only for p = 2; dev/plan.md C253)
  this->constant = lp_unit_ball_log_inverse_volume(n, this->p);
}

arma::mat LpUniformABCKernelFactor::likelihood_evaluate_gradient(const std::string &variable,
                                                                 const Parameters &input) const
{
  return arma::zeros(input[variable].n_elem,1); // (the kernel is flat where it is not 0; it was uninitialised)
}


// (dev/plan.md 3.21) the kernel at distance 0: its largest value
double LpUniformABCKernelFactor::log_maximum(const Parameters &input) const
{
  double scale_constant = (this->scale_variable != "") ? arma::sum(log(arma::colvec(input[this->scale_variable]))) : 0.0;
  return this->evaluate_kernel_given_distance_at(input[this->epsilon_variable][0], 0.0, scale_constant);
}
}
