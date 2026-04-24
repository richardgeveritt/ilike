#ifndef FACTORS_H
#define FACTORS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "distributions.h"
#include "parameters.h"

namespace ilike
{
  /**
   * @file factors.h
   * @brief Defines the Index class.
   *
   * Provides index functionality.
   *
   * @namespace ilike
   * @class Index
   * @brief The index class.
   */


class Index;
class FactorVariables;

class Factors
{
  
public:
  
  /**
   * @brief Performs the factors operation.
   */
  Factors();
  /**
   * @brief Performs the ~factors operation.
   */
  virtual ~Factors();
  
  /**
   * @brief Performs the factors operation.
   *
   * @param another The Index instance to copy from.
   */
  Factors(const Factors &another);
  
  /**
   * @brief Assignment operator for Index.
   *
   * @param another The Index instance to copy from.
   */
  void operator=(const Factors &another);
  /**
   * @brief Creates a deep copy of this Index object.
   *
   * @return The result.
   */
  virtual Factors* duplicate() const=0;
  
  /**
   * @brief Sets the data.
   *
   * @param index The index.
   */
  void set_data(size_t index);
  /**
   * @brief Sets the data.
   *
   * @param index The index.
   */
  virtual void set_data(const Index* index)=0;
  
  /**
   * @brief Performs the change data operation.
   *
   * @param new_data The new data.
   */
  void change_data(Data* new_data);
  
  /**
   * @brief Returns the current data.
   *
   * @return The result.
   */
  virtual Data* get_current_data()=0;
  
  /**
   * @brief Simulates factor variables.
   *
   * @param simulated_parameters The simulated parameters.
   *
   * @return The result.
   */
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
  /**
   * @brief Performs setup given the supplied parameters.
   *
   * @param conditioned_on_parameters The conditioned on parameters.
   */
  virtual void setup(const Parameters &conditioned_on_parameters)=0;
  
protected:
  
  /**
   * @brief Class-specific implementation for change data.
   *
   * @param new_data The new data.
   */
  virtual void specific_change_data(Data* new_data)=0;
  
  /**
   * @brief Copies the state of another Index into this object.
   *
   * @param another The Index instance to copy from.
   */
  void make_copy(const Factors &another);
  
};
}

#endif
