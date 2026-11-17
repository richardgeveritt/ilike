#ifndef ABCKERNELFACTOR_H
#define ABCKERNELFACTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "likelihood_factor.h"
#include "packing_instructions_map.h"

namespace ilike
{

  class ABCKernelFactor : public LikelihoodFactor
  {

  public:
    ABCKernelFactor();

    ABCKernelFactor(const std::vector<std::string> &data_variables_in,
                    const std::string &epsilon_variable_in,
                    Data *data_in);

    ABCKernelFactor(const std::vector<std::string> &data_variables_in,
                    const std::string &epsilon_variable_in,
                    const std::string &scale_variable_in,
                    Data *data_in);

    virtual ~ABCKernelFactor();

    ABCKernelFactor(const ABCKernelFactor &another);

    ABCKernelFactor &operator=(const ABCKernelFactor &another);

    virtual ABCKernelFactor *abc_kernel_factor_duplicate() const = 0;

    virtual void find_distance(const Parameters &input,
                               double &distance,
                               double &scale_constant) const = 0;

    virtual double evaluate_kernel_given_distance(const Parameters &input,
                                                  double distance,
                                                  double scale_constant) const = 0;

    // The same at a given tolerance (the value of input[epsilon_variable][0]),
    // so that many distances can be evaluated with one lookup of the
    // tolerance (dev/plan.md 2.6b). evaluate_kernel_given_distance is this at
    // get_epsilon(input): one copy of the expression.
    virtual double evaluate_kernel_given_distance_at(double epsilon,
                                                     double distance,
                                                     double scale_constant) const = 0;
    double get_epsilon(const Parameters &input) const { return input[this->epsilon_variable][0]; }

  protected:
    void make_copy(const ABCKernelFactor &another);

    std::vector<std::string> data_variables;

    std::string epsilon_variable;

    std::string scale_variable;

    arma::colvec data_colvec;

    PackingInstructionsMap packing_instructions;
  };
}

#endif
