/**
 * @file parameters.h
 * @brief This file contains the declaration and implementation of the Parameters class.
 *
 * The Parameters class is designed to manage and manipulate parameters, which can be either
 * Armadillo matrices or any other type using Boost's any type. The class supports various
 * operations such as deep copying, merging, and accessing parameters by name.
 *
 * Dependencies:
 * - RcppArmadillo
 * - Boost (unordered_map and any)
 *
 * Classes:
 * - Parameters: Manages parameters stored as Armadillo matrices or Boost any types.
 *
 * Typedefs:
 * - vector_parameter_iterator: Iterator for vector parameters.
 * - vector_parameter_const_iterator: Const iterator for vector parameters.
 * - any_parameter_iterator: Iterator for any parameters.
 * - any_parameter_const_iterator: Const iterator for any parameters.
 *
 * Macros:
 * - RCPP_EXPOSED_CLASS(Parameters): Exposes the Parameters class to Rcpp.
 *
 * Public Methods:
 * - Constructors: Default, copy, move, and parameterized constructors.
 * - Destructor.
 * - Operator Overloads: [], (), =, <<, and arithmetic operators.
 * - Accessors: get_colvec, get_rowvec, get_matrices, row, col, rows, cols, min_n_rows, min_n_cols.
 * - Modifiers: merge, merge_with_fixed, deep_copy, deep_copy_nonfixed, self_deep_copy_nonfixed, deep_overwrite_with_variables_in_argument.
 * - Iterators: vector_begin, vector_end, any_begin, any_end.
 * - Size Methods: vector_size, any_size.
 * - Variable Methods: get_vector_variables, get_any_variables, get_nonfixed_vector_variables, get_nonfixed_any_variables, get_variable_n_elems.
 *
 * Protected Members:
 * - vector_parameters: Stores parameters as Armadillo matrices.
 * - any_parameters: Stores parameters as Boost any types.
 *
 * Friend Functions:
 * - operator<<: Outputs the Parameters object to an ostream.
 *
 * Inline Methods:
 * - Constructors, Destructor, Copy/Move Constructors and Assignment Operators.
 * - Operator Overloads: [], (), =, <<, and arithmetic operators.
 * - Accessors and Modifiers.
 * - Iterators and Size Methods.
 * - Variable Methods.
 */

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <unordered_map>
#include <memory>
#include <RcppArmadillo.h>
#include <iostream>
#include <limits>
#include <vector>
#include <string>
#include <memory>
#include <boost/unordered_map.hpp>
#include <boost/any.hpp>

class Parameters;

// A stacked variable (dev/plan.md, C203) is the column vector of several
// variables' elements, named by their names joined by ';' ("a;b;c"); a name
// cannot contain ';'. Proposals over several named variables use it as one
// block.
inline bool is_stacked_name(const std::string &name)
{
  return name.find(';') != std::string::npos;
}

inline std::vector<std::string> split_stacked_name(const std::string &name)
{
  std::vector<std::string> names;
  size_t start = 0;
  while (true)
  {
    size_t end = name.find(';', start);
    names.push_back(name.substr(start, (end == std::string::npos) ? std::string::npos : end - start));
    if (end == std::string::npos)
      break;
    start = end + 1;
  }
  return names;
}

// how a block is named in messages: theta, or [a, b, c] for a stack
inline std::string stacked_display_name(const std::string &name)
{
  if (!is_stacked_name(name))
    return name;
  std::string shown;
  for (const std::string &part : split_stacked_name(name))
    shown = shown + (shown.empty() ? "" : ", ") + part;
  return "[" + shown + "]";
}

// a proposal's dimension against the block it moves (C205): `what` names the
// proposal's parameter, e.g. "covariance" (k x k) or "half-width(s)" (k)
inline void check_proposal_size(const std::string &name, size_t elements, size_t proposal_size,
                                const std::string &what, bool square = true)
{
  if (elements == proposal_size)
    return;
  std::string shown = stacked_display_name(name);
  std::string size = std::to_string(proposal_size);
  if (square)
    Rcpp::stop("The proposal's " + what + " for " + shown + " is " + size + " x " + size +
               ", but " + shown + " has " + std::to_string(elements) + " element(s).");
  Rcpp::stop("The proposal has " + size + " " + what + " for " + shown +
             ", but " + shown + " has " + std::to_string(elements) + " element(s).");
}

// A small map from variable names to values, kept as a vector of
// (name, value) pairs searched linearly (dev/plan.md, Stage 3 (b)). A model
// has a handful of variables, so this is faster than a hash map: no hashing,
// one allocation for the whole map instead of one node per variable, and
// copies are a single vector copy. It offers the operations Parameters uses,
// with std::unordered_map's semantics: insert never overwrites an existing
// name; operator[] inserts a default value for a new name. Iteration is in
// insertion order. Values are shared_ptrs in Parameters, so references to
// the matrices stay valid when the vector grows.
// With more than index_threshold names, a hash index from name to position
// makes lookups O(1) (dev/plan.md, "Many named variables": a term reading D
// named variables would otherwise cost D^2 comparisons). The index is
// shared by copies (which keep the same order of names) and is not changed
// while shared (an insertion of a new name then makes a new one); lookups
// never change it, so concurrent lookups are safe.
// The index lives in a base class, empty for maps that never have many names
// (Indexed = false: the non-numeric variables), so those carry no pointer
// (C211: the index's pointer, copied and destroyed with every Parameters,
// cost small models a few percent of an MCMC iteration).
template <bool Indexed>
struct FlatMapIndex
{
  // (shared by copies: changed only when this map is its only owner)
  std::shared_ptr<std::unordered_map<std::string, size_t>> index;

  bool index_find(const std::string &key, size_t &position) const
  {
    if (this->index == NULL)
      return false;
    auto found = this->index->find(key);
    position = (found == this->index->end()) ? size_t(-1) : found->second;
    return true;
  }
  void index_clear() { this->index.reset(); }

  // after a name is appended: the index, once there are enough names (a new
  // table, since copies may share the old one)
  template <typename Data>
  void index_added(const Data &data, size_t threshold)
  {
    if (data.size() <= threshold)
      return;
    if ((this->index != NULL) && (this->index.use_count() == 1))
    {
      // (this map alone holds the table: extended in place)
      this->index->emplace(data.back().first, data.size() - 1);
      return;
    }
    std::shared_ptr<std::unordered_map<std::string, size_t>> table;
    if (this->index != NULL)
      table = std::make_shared<std::unordered_map<std::string, size_t>>(*this->index);
    else
    {
      table = std::make_shared<std::unordered_map<std::string, size_t>>();
      table->reserve(2*data.size());
      for (size_t k = 0; k + 1 < data.size(); ++k)
        table->emplace(data[k].first, k);
    }
    table->emplace(data.back().first, data.size() - 1);
    this->index = table;
  }
};

template <>
struct FlatMapIndex<false>
{
  bool index_find(const std::string &key, size_t &position) const { return false; }
  void index_clear() {}
  template <typename Data>
  void index_added(const Data &data, size_t threshold) {}
};

template <typename V, bool Indexed = true>
class FlatMap : private FlatMapIndex<Indexed>
{
public:
  typedef std::pair<std::string, V> value_type;
  typedef typename std::vector<value_type>::iterator iterator;
  typedef typename std::vector<value_type>::const_iterator const_iterator;

  iterator begin() { return this->data.begin(); }
  iterator end() { return this->data.end(); }
  const_iterator begin() const { return this->data.begin(); }
  const_iterator end() const { return this->data.end(); }
  size_t size() const { return this->data.size(); }
  bool empty() const { return this->data.empty(); }
  void clear() { this->data.clear(); this->index_clear(); }
  void reserve(size_t n) { this->data.reserve(n); }

  iterator find(const std::string &key)
  {
    size_t position;
    if (this->index_find(key, position))
      return (position == size_t(-1)) ? this->data.end() : this->data.begin() + position;
    for (iterator i = this->data.begin(); i != this->data.end(); ++i)
      if (i->first == key)
        return i;
    return this->data.end();
  }

  const_iterator find(const std::string &key) const
  {
    size_t position;
    if (this->index_find(key, position))
      return (position == size_t(-1)) ? this->data.end() : this->data.begin() + position;
    for (const_iterator i = this->data.begin(); i != this->data.end(); ++i)
      if (i->first == key)
        return i;
    return this->data.end();
  }

  std::pair<iterator, bool> insert(const value_type &value)
  {
    iterator found = this->find(value.first);
    if (found != this->data.end())
      return std::pair<iterator, bool>(found, false);
    this->data.push_back(value);
    this->index_added(this->data, index_threshold);
    return std::pair<iterator, bool>(this->data.end() - 1, true);
  }

  std::pair<iterator, bool> insert(value_type &&value)
  {
    iterator found = this->find(value.first);
    if (found != this->data.end())
      return std::pair<iterator, bool>(found, false);
    this->data.push_back(std::move(value));
    this->index_added(this->data, index_threshold);
    return std::pair<iterator, bool>(this->data.end() - 1, true);
  }

  template <typename InputIterator>
  void insert(InputIterator first, InputIterator last)
  {
    for (; first != last; ++first)
      this->insert(*first);
  }

  V& operator[](const std::string &key)
  {
    iterator found = this->find(key);
    if (found != this->data.end())
      return found->second;
    this->data.emplace_back(key, V());
    this->index_added(this->data, index_threshold);
    return this->data.back().second;
  }

private:
  static const size_t index_threshold = 16;

  std::vector<value_type> data;
};

typedef FlatMap<std::pair<std::shared_ptr<arma::mat>, bool>> VectorParameterMap;
typedef FlatMap<std::pair<std::shared_ptr<boost::any>, bool>, false> AnyParameterMap;

typedef VectorParameterMap::iterator vector_parameter_iterator;
typedef VectorParameterMap::const_iterator vector_parameter_const_iterator;
typedef AnyParameterMap::iterator any_parameter_iterator;
typedef AnyParameterMap::const_iterator any_parameter_const_iterator;

RCPP_EXPOSED_CLASS(Parameters)

using boost::any_cast;

/**
 * @class Parameters
 * @brief A class to manage parameters with various types and functionalities.
 *
 * The Parameters class provides a flexible way to handle parameters, supporting
 * operations such as deep copying, merging, and accessing parameters by name.
 * It supports both Armadillo matrices and Boost any types.
 *
 * @note This class uses Armadillo for matrix operations and Boost for any type handling.
 *
 * @file parameters.h
 */

class Parameters
{

public:
  /**
   * @brief Default constructor.
   */
  Parameters();

  /**
   * @brief Constructor with a variable name and a double value.
   * @param variable_in The name of the variable.
   * @param value_in The double value of the variable.
   */
  Parameters(const std::string &variable_in,
             double value_in);

  /**
   * @brief Constructor with a variable name and an Armadillo matrix.
   * @param variable_in The name of the variable.
   * @param value_in The Armadillo matrix value of the variable.
   */
  Parameters(const std::string &variable_in,
             const arma::mat &value_in);

  /**
   * @brief Destructor.
   */
  ~Parameters();

  /**
   * @brief Accessor operator for non-const objects.
   * @param variable The name of the variable.
   * @return Reference to the Armadillo matrix associated with the variable.
   */
  arma::mat &operator[](const std::string &variable);

  /**
   * @brief Accessor operator for const objects.
   * @param variable The name of the variable.
   * @return The Armadillo matrix associated with the variable (a const
   * reference: bind it to a reference, `const arma::mat &x = p["x"];`, to use
   * it without copying; `arma::mat x = p["x"];` copies as before).
   */
  const arma::mat& operator[](const std::string &variable) const;

  // The matrix of a variable, or NULL if there is none (no error; uses the
  // index of a map with many names, dev/plan.md "Many named variables").
  const arma::mat* find_vector(const std::string &variable) const
  {
    auto found = this->vector_parameters.find(variable);
    return (found != this->vector_parameters.end()) ? found->second.first.get() : NULL;
  }

  /**
   * @brief Accessor operator for a vector of variable names.
   * @param variables The vector of variable names.
   * @return Armadillo column vector containing the values of the specified variables.
   */
  arma::colvec operator[](const std::vector<std::string> &variables) const;

  /**
   * @brief Accessor operator for non-const objects using Boost any type.
   * @param variable The name of the variable.
   * @return Reference to the Boost any type associated with the variable.
   */
  boost::any &operator()(const std::string &variable);

  /**
   * @brief Accessor operator for const objects using Boost any type.
   * @param variable The name of the variable.
   * @return Copy of the Boost any type associated with the variable.
   */
  boost::any operator()(const std::string &variable) const;

  /**
   * @brief Copy constructor.
   * @param another The Parameters object to copy from.
   */
  Parameters(const Parameters &another);

  /**
   * @brief Copy assignment operator.
   * @param another The Parameters object to copy from.
   * @return Reference to the current object.
   */
  Parameters &operator=(const Parameters &another);

  /**
   * @brief Creates a duplicate of the current object.
   * @return Pointer to the duplicated Parameters object.
   */
  Parameters *duplicate() const;

  /**
   * @brief Move constructor.
   * @param another The Parameters object to move from.
   */
  Parameters(Parameters &&another);

  /**
   * @brief Move assignment operator.
   * @param another The Parameters object to move from.
   * @return Reference to the current object.
   */
  Parameters &operator=(Parameters &&another);

  /**
   * @brief Checks if the Parameters object is empty.
   * @return True if the object is empty, false otherwise.
   */
  bool is_empty() const;

  /**
   * @brief Gets an Armadillo column vector for a specified variable.
   * @param variable The name of the variable.
   * @return Armadillo column vector associated with the variable.
   */
  arma::colvec get_colvec(const std::string &variable) const;

  // Writes a proposed value of `variable`, with the shape it has in `shapes`
  // (a column if it is not there); for a stacked name, each variable's part
  // of `value`.
  void set_stacked(const std::string &variable,
                   const arma::colvec &value,
                   const Parameters &shapes);

  /**
   * @brief Gets an Armadillo column vector for a vector of variables.
   * @param variables The vector of variable names.
   * @return Armadillo column vector containing the values of the specified variables.
   */
  arma::colvec get_colvec(const std::vector<std::string> &variables) const;

  /**
   * @brief Gets an Armadillo row vector for a specified variable.
   * @param variable The name of the variable.
   * @return Armadillo row vector associated with the variable.
   */
  arma::rowvec get_rowvec(const std::string &variable) const;

  /**
   * @brief Gets an Armadillo row vector for a vector of variables.
   * @param variables The vector of variable names.
   * @return Armadillo row vector containing the values of the specified variables.
   */
  arma::rowvec get_rowvec(const std::vector<std::string> &variables) const;

  /**
   * @brief Gets a vector of Armadillo matrices for a vector of variables.
   * @param variables The vector of variable names.
   * @return Vector of Armadillo matrices associated with the specified variables.
   */
  std::vector<arma::mat> get_matrices(const std::vector<std::string> &variables) const;

  /**
   * @brief Gets a Parameters object representing a specific row.
   * @param index The index of the row.
   * @return Parameters object representing the specified row.
   */
  Parameters row(size_t index) const;

  /**
   * @brief Gets a Parameters object representing a specific column.
   * @param index The index of the column.
   * @return Parameters object representing the specified column.
   */
  Parameters col(size_t index) const;

  /**
   * @brief Gets a Parameters object representing specific rows.
   * @param indices The indices of the rows.
   * @return Parameters object representing the specified rows.
   */
  Parameters rows(const arma::uvec &indices) const;

  /**
   * @brief Gets a Parameters object representing specific columns.
   * @param indices The indices of the columns.
   * @return Parameters object representing the specified columns.
   */
  Parameters cols(const arma::uvec &indices) const;

  /**
   * @brief Gets the minimum number of rows among all parameters.
   * @return Minimum number of rows.
   */
  size_t min_n_rows() const;

  /**
   * @brief Gets the minimum number of columns among all parameters.
   * @return Minimum number of columns.
   */
  size_t min_n_cols() const;

  /**
   * @brief Merges the current Parameters object with another.
   * @param another The Parameters object to merge with.
   * @return Merged Parameters object.
   */
  Parameters merge(const Parameters &another) const;

  /**
   * @brief Merges the current Parameters object with another, where the argument contains fixed parameters.
   * @param conditioned_on_parameters The Parameters object with fixed parameters.
   */
  void merge_with_fixed(const Parameters &conditioned_on_parameters);

  /**
   * @brief merge_with_fixed with two sets of fixed parameters in turn, with
   * room for both made once.
   */
  void merge_with_fixed(const Parameters &first_conditioned_on_parameters,
                        const Parameters &second_conditioned_on_parameters);

  /**
   * @brief Creates a deep copy of the current Parameters object.
   * @return Deep copied Parameters object.
   */
  Parameters deep_copy() const;

  /**
   * @brief Creates a deep copy of the current Parameters object, with shallow copy of fixed parameters.
   * @return Deep copied Parameters object, with shallow copy of fixed parameters.
   */
  Parameters deep_copy_nonfixed() const;

  /**
   * @brief Performs a deep copy of the non-fixed parameters.
   */
  void self_deep_copy_nonfixed();

  /**
   * @brief Deeply overwrites the current object with variables from another Parameters object.
   * @param new_parameters The Parameters object with new variables.
   */
  void deep_overwrite_with_variables_in_argument(const Parameters &new_parameters);

  /**
   * @brief (C279) Sets one variable as deep_overwrite_with_variables_in_argument
   * does: in place when this object alone owns the matrix and the shape is
   * the same, otherwise in a new matrix (a copy of this Parameters sharing
   * the old one does not see the change).
   */
  void overwrite_variable(const std::string &variable, const arma::mat &value);

  /**
   * @brief Gets an iterator to the beginning of vector parameters.
   * @return Iterator to the beginning of vector parameters.
   */
  vector_parameter_iterator vector_begin();

  /**
   * @brief Gets an iterator to the end of vector parameters.
   * @return Iterator to the end of vector parameters.
   */
  vector_parameter_iterator vector_end();

  /**
   * @brief Gets an iterator to the beginning of any parameters.
   * @return Iterator to the beginning of any parameters.
   */
  any_parameter_iterator any_begin();

  /**
   * @brief Gets an iterator to the end of any parameters.
   * @return Iterator to the end of any parameters.
   */
  any_parameter_iterator any_end();

  /**
   * @brief Gets a const iterator to the beginning of vector parameters.
   * @return Const iterator to the beginning of vector parameters.
   */
  vector_parameter_const_iterator vector_begin() const;

  /**
   * @brief Gets a const iterator to the end of vector parameters.
   * @return Const iterator to the end of vector parameters.
   */
  vector_parameter_const_iterator vector_end() const;

  /**
   * @brief Gets a const iterator to the beginning of any parameters.
   * @return Const iterator to the beginning of any parameters.
   */
  any_parameter_const_iterator any_begin() const;

  /**
   * @brief Gets a const iterator to the end of any parameters.
   * @return Const iterator to the end of any parameters.
   */
  any_parameter_const_iterator any_end() const;

  /**
   * @brief Gets the size of vector parameters.
   * @return Size of vector parameters.
   */
  size_t vector_size() const;

  /**
   * @brief Gets the size of any parameters.
   * @return Size of any parameters.
   */
  size_t any_size() const;

  /**
   * @brief Gets the names of vector variables.
   * @return Vector of names of vector variables.
   */
  std::vector<std::string> get_vector_variables() const;

  /**
   * @brief Gets the names of any variables.
   * @return Vector of names of any variables.
   */
  std::vector<std::string> get_any_variables() const;

  /**
   * @brief Gets the names of non-fixed vector variables.
   * @return Vector of names of non-fixed vector variables.
   */
  std::vector<std::string> get_nonfixed_vector_variables() const;

  /**
   * @brief Gets the names of non-fixed any variables.
   * @return Vector of names of non-fixed any variables.
   */
  std::vector<std::string> get_nonfixed_any_variables() const;

  /**
   * @brief Gets the number of elements for each variable in a vector of variables.
   * @param variables The vector of variable names.
   * @return Vector of sizes corresponding to each variable.
   */
  std::vector<size_t> get_variable_n_elems(const std::vector<std::string> &variables) const;

  /**
   * @brief Output stream operator for Parameters.
   * @param os The output stream.
   * @param p The Parameters object.
   * @return Reference to the output stream.
   */
  friend std::ostream &operator<<(std::ostream &os, const Parameters &p);

protected:
  /**
   * @brief A map that associates parameter names with their corresponding matrix and a boolean flag.
   *
   * This unordered map uses a string as the key, which represents the name of the parameter.
   * The value is a pair consisting of:
   * - A shared pointer to an Armadillo matrix (arma::mat), which holds the parameter data.
   * - A boolean flag indicating whether the parameter is treated as "fixed". When Parameters is copied, the fixed parameters only copy the shared pointer, while the non-fixed parameters perform a deep copy.
   */
  VectorParameterMap vector_parameters;

  /**
   * @brief A map that associates parameter names with their values and a flag indicating if they are set.
   *
   * This unordered map uses a string as the key to represent the parameter name.
   * The value is a pair consisting of:
   * - A shared pointer to a boost::any object, which can hold any type of parameter value.
   * - A boolean flag indicating whether the parameter is treated as "fixed". When Parameters is copied, the fixed parameters only copy the shared pointer, while the non-fixed parameters perform a deep copy.
   */
  AnyParameterMap any_parameters;

  /**
   * @brief Copies the content of another Parameters object.
   * @param another The Parameters object to copy from.
   */
  void make_copy(const Parameters &another);

  /**
   * @brief Moves the content of another Parameters object.
   * @param another The Parameters object to move from.
   */
  void make_copy(Parameters &&another);
};

typedef Parameters Data;
typedef Parameters MatrixList;

#include <RcppArmadillo.h>

inline Parameters::Parameters()
{
  this->vector_parameters.clear();
  this->any_parameters.clear();
}

inline Parameters::Parameters(const std::string &variable_in,
                              double value_in)
{
  (*this)[variable_in] = value_in;
}

inline Parameters::Parameters(const std::string &variable_in,
                              const arma::mat &value_in)
{
  (*this)[variable_in] = value_in;
}

inline Parameters::~Parameters()
{
}

// Copy constructor for the Parameters class.
inline Parameters::Parameters(const Parameters &another)
{
  this->make_copy(another);
}

inline Parameters &Parameters::operator=(const Parameters &another)
{
  if (this == &another)
  { // if a==a
    return *this;
  }

  this->make_copy(another);
  return *this;
}

// Move constructor for the Parameters class.
inline Parameters::Parameters(Parameters &&another)
{
  this->make_copy(std::move(another));
}

inline Parameters &Parameters::operator=(Parameters &&another)
{
  if (this == &another)
  { // if a==a
    return *this;
  }

  this->make_copy(std::move(another));
  return *this;
}

inline Parameters *Parameters::duplicate() const
{
  return (new Parameters(*this));
}

inline bool Parameters::is_empty() const
{
  if ((this->vector_parameters.size() == 0) && (this->any_parameters.size() == 0))
    return TRUE;
  else
    return FALSE;
}

inline arma::colvec Parameters::get_colvec(const std::string &variable) const
{
  const arma::mat* found = this->find_vector(variable);
  if (found != NULL)
    return arma::vectorise(*found);
  // a stacked name is looked at only when the name is not a variable
  if (is_stacked_name(variable))
    return this->get_colvec(split_stacked_name(variable));
  return arma::vectorise((*this)[variable]); // the "not found" error
}

inline void Parameters::set_stacked(const std::string &variable,
                                    const arma::colvec &value,
                                    const Parameters &shapes)
{
  if (!is_stacked_name(variable))
  {
    // the shape the variable has at the current point (C205), else a column
    const arma::mat* shape = shapes.find_vector(variable);
    if ((shape != NULL) && (shape->n_elem == value.n_elem) && (shape->n_cols != 1))
      (*this)[variable] = arma::reshape(value, shape->n_rows, shape->n_cols);
    else
      (*this)[variable] = value;
    return;
  }
  std::vector<std::string> names = split_stacked_name(variable);
  arma::uword offset = 0;
  for (const std::string &name : names)
  {
    const arma::mat &shape = shapes[name];
    if (shape.n_elem == 0)
    {
      (*this)[name] = arma::mat(shape.n_rows, shape.n_cols);
      continue;
    }
    if (offset + shape.n_elem > value.n_elem)
      Rcpp::stop("Parameters::set_stacked: the proposed value of " + variable + " has " +
                 std::to_string(value.n_elem) + " elements, fewer than its variables.");
    (*this)[name] = arma::reshape(value.subvec(offset, offset + shape.n_elem - 1), shape.n_rows, shape.n_cols);
    offset = offset + shape.n_elem;
  }
  if (offset != value.n_elem)
    Rcpp::stop("Parameters::set_stacked: the proposed value of " + variable + " has " +
               std::to_string(value.n_elem) + " elements, but its variables have " + std::to_string(offset) + ".");
}

inline arma::colvec Parameters::get_colvec(const std::vector<std::string> &variables) const
{
  arma::colvec concatenated_vector;
  for (auto i = variables.begin();
       i != variables.end();
       ++i)
  {
    if (i == variables.begin())
      concatenated_vector = this->get_colvec(*i);
    else
      concatenated_vector = join_cols(concatenated_vector, this->get_colvec(*i));
  }
  return concatenated_vector;
}

// column by column, as get_colvec (and as R stores a matrix), so outputs,
// population covariances and proposals agree on the order of a matrix's
// elements
inline arma::rowvec Parameters::get_rowvec(const std::string &variable) const
{
  // (vectorise(m, 1) would concatenate the rows: row by row)
  return arma::vectorise((*this)[variable]).t();
}

inline arma::rowvec Parameters::get_rowvec(const std::vector<std::string> &variables) const
{
  arma::rowvec concatenated_vector;
  for (auto i = variables.begin();
       i != variables.end();
       ++i)
  {
    if (i == variables.begin())
      concatenated_vector = this->get_rowvec(*i);
    else
      concatenated_vector = join_rows(concatenated_vector, this->get_rowvec(*i));
  }
  return concatenated_vector;
}

inline std::vector<arma::mat> Parameters::get_matrices(const std::vector<std::string> &variables) const
{
  std::vector<arma::mat> output;
  output.reserve(variables.size());
  for (auto i = variables.begin();
       i != variables.end();
       ++i)
  {
    output.push_back((*this)[*i]);
  }
  return output;
}

inline Parameters Parameters::row(size_t index) const
{
  Parameters output;
  for (auto i = this->vector_begin(); i != this->vector_end(); ++i)
  {
    output.vector_parameters.insert({i->first, std::pair<std::shared_ptr<arma::mat>, bool>(std::make_shared<arma::mat>(i->second.first->row(index)), i->second.second)});
  }
  return output;
}

inline Parameters Parameters::col(size_t index) const
{
  Parameters output;
  for (auto i = this->vector_begin(); i != this->vector_end(); ++i)
  {
    output.vector_parameters.insert({i->first, std::pair<std::shared_ptr<arma::mat>, bool>(std::make_shared<arma::mat>(i->second.first->col(index)), i->second.second)});
  }
  return output;
}

inline Parameters Parameters::rows(const arma::uvec &indices) const
{
  Parameters output;
  for (auto i = this->vector_begin(); i != this->vector_end(); ++i)
  {
    output.vector_parameters.insert({i->first, std::pair<std::shared_ptr<arma::mat>, bool>(std::make_shared<arma::mat>(i->second.first->rows(indices)), i->second.second)});
  }
  return output;
}

inline Parameters Parameters::cols(const arma::uvec &indices) const
{
  Parameters output;
  for (auto i = this->vector_begin(); i != this->vector_end(); ++i)
  {
    output.vector_parameters.insert({i->first, std::pair<std::shared_ptr<arma::mat>, bool>(std::make_shared<arma::mat>(i->second.first->cols(indices)), i->second.second)});
  }
  return output;
}

inline size_t Parameters::min_n_rows() const
{
  size_t min = std::numeric_limits<size_t>::max();
  for (auto i = this->vector_begin(); i != this->vector_end(); ++i)
  {
    if (i->second.first->n_rows < min)
      min = i->second.first->n_rows;
  }
  return min;
}

inline size_t Parameters::min_n_cols() const
{
  size_t min = std::numeric_limits<size_t>::max();
  for (auto i = this->vector_begin(); i != this->vector_end(); ++i)
  {
    if (i->second.first->n_rows < min)
      min = i->second.first->n_cols;
  }
  return min;
}

inline void Parameters::make_copy(const Parameters &another)
{
  this->vector_parameters = another.vector_parameters;
  this->any_parameters = another.any_parameters;
}

inline void Parameters::make_copy(Parameters &&another)
{
  this->vector_parameters = std::move(another.vector_parameters);
  this->any_parameters = std::move(another.any_parameters);

  // leave the source empty without constructing temporary maps
  another.vector_parameters.clear();
  another.any_parameters.clear();
}

inline arma::mat &Parameters::operator[](const std::string &variable)
{
  auto found = this->vector_parameters.find(variable);

  // if doesn't yet exist, make shared pointer, and set to be not fixed
  if (found == this->vector_parameters.end())
  {
    return *this->vector_parameters.insert({variable, std::pair<std::shared_ptr<arma::mat>, bool>(std::make_shared<arma::mat>(), false)}).first->second.first;
  }
  else
  {
    if (found->second.second == false)
    {
      // if does exist, and is not fixed, just pass reference to memory that the shard pointer points to
      return *found->second.first;
    }
    else
    {
      // if does exist, and is fixed, throw error
      Rcpp::stop("Parameters::operator[] - parameter is fixed and cannot be changed.");
    }
  }
}

inline const arma::mat& Parameters::operator[](const std::string &variable) const
{
  auto found = this->vector_parameters.find(variable);
  if (found != this->vector_parameters.end())
    return (*found->second.first);
  else
    Rcpp::stop("Parameters::operator[]: variable '" + variable + "' not found in Parameters.");
}

inline arma::colvec Parameters::operator[](const std::vector<std::string> &variables) const
{
  arma::mat concatenated_matrix;
  for (std::vector<std::string>::const_iterator i = variables.begin();
       i != variables.end();
       ++i)
  {
    auto found = this->vector_parameters.find(*i);
    if (i == variables.begin())
    {
      if (found != this->vector_parameters.end())
        concatenated_matrix = arma::vectorise(*found->second.first);
      else
        Rcpp::stop("Parameters::operator[]: variable not found in Parameters.");
    }
    else
    {
      if (found != this->vector_parameters.end())
        concatenated_matrix = join_cols(concatenated_matrix, arma::vectorise(*found->second.first));
      else
        Rcpp::stop("Parameters::operator[]: variable not found in Parameters.");
    }
  }

  return concatenated_matrix;
}

inline boost::any &Parameters::operator()(const std::string &variable)
{
  auto found = this->any_parameters.find(variable);

  // if doesn't yet exist, make shared pointer, and set to be not fixed
  if (found == this->any_parameters.end())
  {
    return *this->any_parameters.insert({variable, std::pair<std::shared_ptr<boost::any>, bool>(std::make_shared<boost::any>(), false)}).first->second.first;
  }
  else
  {
    if (found->second.second == false)
    {
      // if does exist, and is not fixed, just pass reference to memory that the shard pointer points to
      return *found->second.first;
    }
    else
    {
      // if does exist, and is fixed, throw error
      Rcpp::stop("Parameters::operator() - parameter is fixed and cannot be changed.");
    }
  }
}

inline boost::any Parameters::operator()(const std::string &variable) const
{
  auto found = this->any_parameters.find(variable);

  if (found != this->any_parameters.end())
    return (*found->second.first);
  else
    Rcpp::stop("Parameters::operator(): variable not found in Parameters.");
}

inline Parameters Parameters::merge(const Parameters &another) const
{
  Parameters output;
  output.vector_parameters = this->vector_parameters;
  output.any_parameters = this->any_parameters;

  output.vector_parameters.insert(another.vector_parameters.begin(), another.vector_parameters.end());

  output.any_parameters.insert(another.any_parameters.begin(), another.any_parameters.end());

  return output;
}

inline void Parameters::merge_with_fixed(const Parameters &conditioned_on_parameters)
{
  // room for the new variables in one allocation (more than needed if some
  // names are already here)
  if (conditioned_on_parameters.vector_parameters.size() > 0)
    this->vector_parameters.reserve(this->vector_parameters.size() + conditioned_on_parameters.vector_parameters.size());
  if (conditioned_on_parameters.any_parameters.size() > 0)
    this->any_parameters.reserve(this->any_parameters.size() + conditioned_on_parameters.any_parameters.size());

  for (auto i = conditioned_on_parameters.vector_begin(); i != conditioned_on_parameters.vector_end(); ++i)
  {
    this->vector_parameters[i->first] = std::pair<std::shared_ptr<arma::mat>, bool>(i->second.first, true);
  }

  for (auto i = conditioned_on_parameters.any_begin(); i != conditioned_on_parameters.any_end(); ++i)
  {
    this->any_parameters[i->first] = std::pair<std::shared_ptr<boost::any>, bool>(i->second.first, true);
  }
}

inline void Parameters::merge_with_fixed(const Parameters &first_conditioned_on_parameters,
                                         const Parameters &second_conditioned_on_parameters)
{
  this->vector_parameters.reserve(this->vector_parameters.size() +
                                  first_conditioned_on_parameters.vector_parameters.size() +
                                  second_conditioned_on_parameters.vector_parameters.size());
  this->any_parameters.reserve(this->any_parameters.size() +
                               first_conditioned_on_parameters.any_parameters.size() +
                               second_conditioned_on_parameters.any_parameters.size());
  this->merge_with_fixed(first_conditioned_on_parameters);
  this->merge_with_fixed(second_conditioned_on_parameters);
}

inline Parameters Parameters::deep_copy() const
{
  Parameters output;

  for (auto i = this->vector_begin();
       i != this->vector_end();
       ++i)
  {
    if (i->second.second == true) // fixed
    {
      output.vector_parameters.insert({i->first, std::pair<std::shared_ptr<arma::mat>, bool>(std::make_shared<arma::mat>(*i->second.first), true)});
    }
    else // non-fixed
    {
      output.vector_parameters.insert({i->first, std::pair<std::shared_ptr<arma::mat>, bool>(std::make_shared<arma::mat>(*i->second.first), false)});
    }
  }

  for (auto i = this->any_begin();
       i != this->any_end();
       ++i)
  {
    if (i->second.second == true) // fixed
    {
      output.any_parameters.insert({i->first, std::pair<std::shared_ptr<boost::any>, bool>(std::make_shared<boost::any>(*i->second.first), true)});
    }
    else // non-fixed
    {
      output.any_parameters.insert({i->first, std::pair<std::shared_ptr<boost::any>, bool>(std::make_shared<boost::any>(*i->second.first), false)});
    }
  }

  return output;
}

inline Parameters Parameters::deep_copy_nonfixed() const
{
  Parameters output;

  for (auto i = this->vector_begin();
       i != this->vector_end();
       ++i)
  {
    if (i->second.second == true) // fixed, shallow copy
    {
      output.vector_parameters.insert({i->first, std::pair<std::shared_ptr<arma::mat>, bool>(i->second.first, true)});
    }
    else // non-fixed, deep copy
    {
      output.vector_parameters.insert({i->first, std::pair<std::shared_ptr<arma::mat>, bool>(std::make_shared<arma::mat>(*i->second.first), false)});
    }
  }

  for (auto i = this->any_begin();
       i != this->any_end();
       ++i)
  {
    if (i->second.second == true) // fixed, shallow copy
    {
      output.any_parameters.insert({i->first, std::pair<std::shared_ptr<boost::any>, bool>(i->second.first, true)});
    }
    else // non-fixed, deep copy
    {
      output.any_parameters.insert({i->first, std::pair<std::shared_ptr<boost::any>, bool>(std::make_shared<boost::any>(*i->second.first), false)});
    }
  }

  return output;
}

inline void Parameters::self_deep_copy_nonfixed()
{

  for (auto i = this->vector_begin();
       i != this->vector_end();
       ++i)
  {
    if (i->second.second == false) // not fixed, deep copy
    {
      i->second.first = std::make_shared<arma::mat>(*i->second.first);
    }
  }

  for (auto i = this->any_begin();
       i != this->any_end();
       ++i)
  {
    if (i->second.second == false) // not fixed, deep copy
    {
      i->second.first = std::make_shared<boost::any>(*i->second.first);
    }
  }
}

inline void Parameters::overwrite_variable(const std::string &variable, const arma::mat &value)
{
  auto found = this->vector_parameters.find(variable);
  if (found == this->vector_parameters.end())
  {
    this->vector_parameters.insert({variable, std::pair<std::shared_ptr<arma::mat>, bool>(std::make_shared<arma::mat>(value), false)});
    return;
  }
  if (found->second.second == true)
    Rcpp::stop("Parameters::overwrite_variable - parameter is fixed and cannot be changed.");
  std::shared_ptr<arma::mat> &old_matrix = found->second.first;
  if ((old_matrix.use_count() == 1) && (old_matrix->n_rows == value.n_rows) && (old_matrix->n_cols == value.n_cols))
    *old_matrix = value;
  else
    old_matrix = std::make_shared<arma::mat>(value);
}

inline void Parameters::deep_overwrite_with_variables_in_argument(const Parameters &new_parameters)
{
  for (auto i = new_parameters.vector_begin(); i != new_parameters.vector_end(); ++i)
  {
    auto found = this->vector_parameters.find(i->first);

    // if doesn't yet exist, make shared pointer, and set to be not fixed
    if (found == this->vector_parameters.end())
    {
      this->vector_parameters.insert({i->first, std::pair<std::shared_ptr<arma::mat>, bool>(std::make_shared<arma::mat>(*i->second.first), false)});
    }
    else
    {
      if (found->second.second == false)
      {
        // If does exist, and is not fixed: replace the values. A new matrix is
        // needed only when the old one is shared with another Parameters
        // object (which must not see the change) or has another shape; a
        // matrix owned by this object alone is overwritten in place, which
        // gives the same values without an allocation.
        std::shared_ptr<arma::mat> &old_matrix = found->second.first;
        const arma::mat &new_matrix = *i->second.first;
        if ( (old_matrix.use_count() == 1) &&
             (old_matrix->n_rows == new_matrix.n_rows) &&
             (old_matrix->n_cols == new_matrix.n_cols) )
          *old_matrix = new_matrix;
        else
          old_matrix = std::make_shared<arma::mat>(new_matrix);
      }
      else
      {
        // if does exist, and is fixed, throw error
        Rcpp::stop("Parameters::deep_overwrite_with_variables_in_argument - parameter is fixed and cannot be changed.");
      }
    }
  }

  for (auto i = new_parameters.any_begin(); i != new_parameters.any_end(); ++i)
  {
    auto found = this->any_parameters.find(i->first);

    // if doesn't yet exist, make shared pointer, and set to be not fixed
    if (found == this->any_parameters.end())
    {
      this->any_parameters.insert({i->first, std::pair<std::shared_ptr<boost::any>, bool>(std::make_shared<boost::any>(*i->second.first), false)});
    }
    else
    {
      if (found->second.second == false)
      {
        // if does exist, and is not fixed, just pass reference to memory that the shard pointer points to
        found->second.first = std::make_shared<boost::any>(*i->second.first);
      }
      else
      {
        // if does exist, and is fixed, throw error
        Rcpp::stop("Parameters::deep_overwrite_with_variables_in_argument - parameter is fixed and cannot be changed.");
      }
    }
  }
}

inline vector_parameter_iterator Parameters::vector_begin()
{
  return this->vector_parameters.begin();
}

inline vector_parameter_iterator Parameters::vector_end()
{
  return this->vector_parameters.end();
}

inline any_parameter_iterator Parameters::any_begin()
{
  return this->any_parameters.begin();
}

inline any_parameter_iterator Parameters::any_end()
{
  return this->any_parameters.end();
}

inline vector_parameter_const_iterator Parameters::vector_begin() const
{
  return this->vector_parameters.begin();
}

inline vector_parameter_const_iterator Parameters::vector_end() const
{
  return this->vector_parameters.end();
}

inline any_parameter_const_iterator Parameters::any_begin() const
{
  return this->any_parameters.begin();
}

inline any_parameter_const_iterator Parameters::any_end() const
{
  return this->any_parameters.end();
}

inline size_t Parameters::vector_size() const
{
  return this->vector_parameters.size();
}

inline size_t Parameters::any_size() const
{
  return this->any_parameters.size();
}

inline std::vector<std::string> Parameters::get_vector_variables() const
{
  std::vector<std::string> variables;
  for (auto it = this->vector_parameters.begin(); it != this->vector_parameters.end(); ++it)
  {
    variables.push_back(it->first);
  }
  return variables;
}

inline std::vector<std::string> Parameters::get_any_variables() const
{
  std::vector<std::string> variables;
  for (auto it = this->any_parameters.begin(); it != this->any_parameters.end(); ++it)
  {
    variables.push_back(it->first);
  }
  return variables;
}

inline std::vector<std::string> Parameters::get_nonfixed_vector_variables() const
{
  std::vector<std::string> variables;
  for (auto it = this->vector_parameters.begin(); it != this->vector_parameters.end(); ++it)
  {
    if (it->second.second == false)
      variables.push_back(it->first);
  }
  return variables;
}

inline std::vector<std::string> Parameters::get_nonfixed_any_variables() const
{
  std::vector<std::string> variables;
  for (auto it = this->any_parameters.begin(); it != this->any_parameters.end(); ++it)
  {
    if (it->second.second == false)
      variables.push_back(it->first);
  }
  return variables;
}

inline std::vector<size_t> Parameters::get_variable_n_elems(const std::vector<std::string> &variables) const
{
  std::vector<size_t> n_elems;
  n_elems.reserve(variables.size());
  for (size_t i = 0; i < variables.size(); ++i)
  {
    n_elems.push_back((*this)[variables[i]].n_elem);
  }
  return n_elems;
}

inline Parameters pow(const Parameters &p, double power)
{
  Parameters output;
  for (auto it = p.vector_begin(); it != p.vector_end(); ++it)
  {
    output[it->first] = arma::pow(p[it->first], power);
  }

  return output;
}

inline Parameters operator*(double scale, const Parameters &p)
{
  Parameters output = p;
  for (auto it = p.vector_begin(); it != p.vector_end(); ++it)
  {
    output[it->first] = scale * p[it->first];
  }
  return output;
}

inline std::ostream &operator<<(std::ostream &os, const Parameters &p)
{
  for (auto it = p.vector_parameters.begin(); it != p.vector_parameters.end(); ++it)
  {
    if (it->first != "")
    {
      os << it->first << "=";
      for (arma::mat::const_iterator i = it->second.first->begin(); i != it->second.first->end(); ++i)
      {
        if (i == it->second.first->begin())
          os << *i;
        else
          os << ", " << *i;
      }
      os << ";";
    }
  }

  for (auto it2 = p.any_parameters.begin(); it2 != p.any_parameters.end(); ++it2)
  {
    if (it2->first != "")
    {
      os << it2->first << "=";
      os << it2->second.first;
      os << std::endl;
    }
  }

  return os;
}

#endif
