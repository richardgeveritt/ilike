#ifndef RECIPE_VALUE_H
#define RECIPE_VALUE_H

#include <RcppArmadillo.h>

#include <initializer_list>
#include <memory>
#include <string>
#include <vector>

namespace ilike
{

// A reference to a function used by a recipe: the language it is written in,
// its name, and (once compiled) a pointer to the place where the C++
// function pointer is stored (for functions compiled by the R binding, the
// address held by an R XPtr). The engine only ever calls the C++ function.
struct FunctionReference
{
  std::string language;
  std::string name;
  const void *pointer = nullptr;
};

// The intermediate representation of a recipe (dev/plan.md 3.14, I1): a tree
// of values, independent of the language that produced it. A value is one
// of: null; a vector of logicals, integers, reals or strings (a scalar is a
// vector of length 1); a real matrix; an ordered list of values, each with a
// name (possibly empty); a function reference; or a value of another kind
// that the core does not interpret (kept so that recipes convert without
// loss, refused if used). Values are cheap to copy (shared nodes); lists are
// built with push_back.
class RecipeValue
{

public:

  enum class Kind { Null, Logical, Integer, Real, String, Matrix, List, Function, Other };

  RecipeValue();
  RecipeValue(double value);
  RecipeValue(int value);
  RecipeValue(size_t value);
  RecipeValue(bool value);
  RecipeValue(const char *value);
  RecipeValue(const std::string &value);
  RecipeValue(const std::vector<double> &values);
  RecipeValue(const std::vector<int> &values);
  RecipeValue(const std::vector<bool> &values);
  RecipeValue(const std::vector<std::string> &values);
  RecipeValue(const arma::mat &value);
  RecipeValue(const arma::colvec &value);
  RecipeValue(const FunctionReference &value);
  RecipeValue(const std::vector<std::vector<double>> &values); // a list of real vectors
  RecipeValue(const RecipeValue &another) = default;
  RecipeValue &operator=(const RecipeValue &another) = default;

  // an empty list, or a list of unnamed values
  static RecipeValue list();
  static RecipeValue list(std::initializer_list<RecipeValue> values);
  static RecipeValue other(const std::string &description);
  static RecipeValue logicals(const std::vector<int> &values);
  static RecipeValue integers(const std::vector<int> &values);

  void push_back(const RecipeValue &value, const std::string &name = "");

  Kind kind() const;
  std::string kind_name() const;
  bool is_null() const;
  bool is_list() const;
  bool is_function() const;
  bool is_string() const;
  bool is_real() const;          // reals (vector or matrix)
  bool is_numeric() const;       // reals, integers or logicals (vector or matrix)
  bool is_matrix() const;
  bool is_vector() const;        // an atomic vector, a matrix or a list (as R's isVector)

  // the length of a vector or list, the number of elements of a matrix
  size_t size() const;
  size_t length() const { return this->size(); }

  // lists
  bool contains(const std::string &name) const;
  const RecipeValue &operator[](const std::string &name) const;
  // (an exact match for string literals: otherwise value["key"] could also
  // be read as the built-in "key"[value], through a numeric conversion)
  const RecipeValue &operator[](const char *name) const { return (*this)[std::string(name)]; }
  const RecipeValue &operator[](size_t index) const;
  const RecipeValue &operator[](int index) const;
  const RecipeValue &operator[](long index) const { return (*this)[int(index)]; }

  // replace the element with this name, or append it
  void set(const std::string &name, const RecipeValue &value);
  std::vector<std::string> names() const;

  // conversions (with R's rules: the first element of a vector; a vector as
  // a one-column matrix; integers and logicals as reals)
  double as_double() const;
  int as_int() const;
  size_t as_size_t() const;
  bool as_bool() const;
  std::string as_string() const;
  arma::mat as_matrix() const;
  arma::colvec as_colvec() const;
  arma::rowvec as_rowvec() const;
  std::vector<double> as_doubles() const;
  std::vector<int> as_ints() const;
  std::vector<size_t> as_size_ts() const;
  std::vector<std::string> as_strings() const;
  const FunctionReference &as_function() const;

  operator double() const { return this->as_double(); }
  operator bool() const { return this->as_bool(); }
  operator int() const { return this->as_int(); }
  operator size_t() const { return this->as_size_t(); }
  operator std::string() const { return this->as_string(); }
  operator arma::mat() const { return this->as_matrix(); }
  operator arma::colvec() const { return this->as_colvec(); }
  operator std::vector<double>() const { return this->as_doubles(); }
  operator std::vector<int>() const { return this->as_ints(); }
  operator std::vector<size_t>() const { return this->as_size_ts(); }
  operator std::vector<std::string>() const { return this->as_strings(); }
  operator std::vector<std::vector<double>>() const { return this->as_list_of_doubles(); }

  // a list of numeric vectors
  std::vector<std::vector<double>> as_list_of_doubles() const;

  // the compiled C++ function a function reference points to
  template <class FunctionPointer>
  FunctionPointer function_pointer() const
  {
    const FunctionReference &reference = this->as_function();
    if (reference.pointer == nullptr)
      Rcpp::stop("Recipe: the function " + reference.name + " has not been compiled.");
    return *static_cast<const FunctionPointer*>(reference.pointer);
  }

  // JSON, for inspection and tests (function pointers are not written)
  std::string to_json() const;

private:

  struct Node
  {
    Kind kind = Kind::Null;
    std::vector<int> logicals;    // also integers
    std::vector<double> reals;
    std::vector<std::string> strings;
    arma::mat matrix;
    std::vector<RecipeValue> items;
    std::vector<std::string> keys;
    FunctionReference function;
    std::string description;      // Other
  };

  std::shared_ptr<Node> node;

  explicit RecipeValue(std::shared_ptr<Node> node_in);
  Node &writable();
  void require(bool condition, const std::string &wanted) const;
  void write_json(std::string &out) const;

};

}

#endif
