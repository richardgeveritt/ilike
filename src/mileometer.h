#ifndef MILEOMETER_H
#define MILEOMETER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

// The point of this is to iterate through a vector of integers in the following way:
// 0 0 0 0
// 0 0 0 1
// 0 0 0 2
// 0 0 1 0
// 0 0 1 1
// 0 0 1 2
// 0 1 0 0
// 0 1 0 1
// etc.

#include <vector>
#include <stdexcept>

namespace ilike
{
  /**
   * @file mileometer.h
   * @brief Defines the Sequencer class.
   *
   * Manages the sequence of tempering or annealing steps in an SMC run. Determines the schedule of bridging distributions.
   *
   * @namespace ilike
   * @class Sequencer
   * @brief The sequencer class.
   */


class Sequencer;
class EnsembleSequencer;

class Mileometer
{
public:
  
  /**
   * @brief Performs the mileometer operation.
   */
  Mileometer();
  /**
   * @brief Performs the ~mileometer operation.
   */
  virtual ~Mileometer();
  
  /**
   * @brief Performs the mileometer operation.
   *
   * @param limitsin The limitsin.
   */
  Mileometer(const std::vector<size_t> &limitsin);
  
  // Everything you need to copy the class.
  /**
   * @brief Performs the mileometer operation.
   *
   * @param another The Sequencer instance to copy from.
   */
  Mileometer(const Mileometer &another);
  /**
   * @brief Creates a deep copy of this Sequencer object.
   *
   * @return The result.
   */
  Mileometer* duplicate() const;
  /**
   * @brief Assignment operator for Sequencer.
   *
   * @param another The Sequencer instance to copy from.
   */
  Mileometer& operator=(const Mileometer &another);
  
  /**
   * @brief Performs the mileometer operation.
   *
   * @param another The Sequencer instance to copy from.
   */
  Mileometer(Mileometer &&another);
  /**
   * @brief Assignment operator for Sequencer.
   *
   * @param another The Sequencer instance to copy from.
   */
  Mileometer& operator=(Mileometer &&another);
  
  // Get the current index.
  /**
   * @brief Returns the current index.
   *
   * @return The result.
   */
  std::vector<size_t> get_current_index() const;
  
  /**
   * @brief Returns the current values.
   *
   * @param values The values.
   *
   * @return The result.
   */
  std::vector<double> get_current_values(const std::vector< std::vector<double> > &values) const;
  
  /**
   * @brief Performs the operator[] operation.
   *
   * @param i The i.
   *
   * @return The result.
   */
  size_t operator[](const size_t &i) const;
  
  /**
   * @brief Performs the size operation.
   *
   * @return The result.
   */
  size_t size() const;
  
  /**
   * @brief Performs the back operation.
   *
   * @return The result.
   */
  size_t back() const;
  
  /**
   * @brief Returns the previous index.
   *
   * @return The result.
   */
  size_t get_previous_index() const;
  /**
   * @brief Returns the next index.
   *
   * @return The result.
   */
  size_t get_next_index() const;
  
  /**
   * @brief Performs the at start operation.
   *
   * @return The result.
   */
  bool at_start() const;
  
  // Increment current.
  /**
   * @brief Performs the increment operation.
   */
  void increment();
  
  /**
   * @brief Performs the reset final dimension operation.
   *
   * @param new_number_of_states The new number of states.
   */
  void reset_final_dimension(size_t new_number_of_states);
  
  /**
   * @brief Performs the reset operation.
   */
  void reset();
  
protected: // Things that can be accessed in this class and subclasses.
  
  friend Sequencer;
  friend EnsembleSequencer;
  
  /** @brief The current index. */
  std::vector<size_t> current_index;
  
  /** @brief The limits. */
  std::vector<size_t> limits;
  
  // Increment current dimension.
  /**
   * @brief Performs the increment operation.
   *
   * @param dimension The dimension.
   */
  void increment(size_t dimension);
  
  /**
   * @brief Copies the state of another Sequencer into this object.
   *
   * @param another The Sequencer instance to copy from.
   */
  void make_copy(const Mileometer &another);
  /**
   * @brief Copies the state of another Sequencer into this object.
   *
   * @param another The Sequencer instance to copy from.
   */
  void make_copy(Mileometer &&another);
  
private: // Things that can be accessed only by this class.
  
};
}

#endif
