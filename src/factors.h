#ifndef FACTORS_H
#define FACTORS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "distributions.h"
#include "parameters.h"

class Index;
class FactorVariables;

class Factors
{

public:

  Factors();
  virtual ~Factors();
  
  Factors(const Factors &another);

  void operator=(const Factors &another);
  virtual Factors* duplicate() const=0;
  
  void set_data(size_t index);
  virtual void set_data(const Index* index)=0;
  
  void change_data(Data* new_data);
  
  virtual Data* get_current_data()=0;
  
  virtual FactorVariables* simulate_factor_variables(const Parameters &simulated_parameters) const=0;
  
  /*
  virtual FactorVariables* simulate_factor_variables(const Parameters &simulated_parameters,
                                                     const Parameters &conditioned_on_parameters)=0;
   
  */
  
  virtual FactorVariables* subsample_simulate_factor_variables(const Parameters &simulated_parameters) const=0;
  
  /*
  virtual FactorVariables* subsample_simulate_factor_variables(const Parameters &simulated_parameters,
                                                               const Parameters &conditioned_on_parameters)=0;
  */
  
  virtual void setup()=0;
  virtual void setup(const Parameters &conditioned_on_parameters)=0;

protected:
  
  virtual void specific_change_data(Data* new_data)=0;

  void make_copy(const Factors &another);

};

#endif
