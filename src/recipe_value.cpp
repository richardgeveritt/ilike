#include <cmath>
#include <cstdio>
#include <limits>
#include "recipe_value.h"

namespace ilike
{

RecipeValue::RecipeValue()
:node(std::make_shared<Node>())
{
}

RecipeValue::RecipeValue(std::shared_ptr<Node> node_in)
:node(node_in)
{
}

RecipeValue::RecipeValue(double value)
:node(std::make_shared<Node>())
{
  this->node->kind = Kind::Real;
  this->node->reals.push_back(value);
}

RecipeValue::RecipeValue(int value)
:node(std::make_shared<Node>())
{
  this->node->kind = Kind::Integer;
  this->node->logicals.push_back(value);
}

RecipeValue::RecipeValue(size_t value)
:node(std::make_shared<Node>())
{
  this->node->kind = Kind::Real;
  this->node->reals.push_back(double(value));
}

RecipeValue::RecipeValue(bool value)
:node(std::make_shared<Node>())
{
  this->node->kind = Kind::Logical;
  this->node->logicals.push_back(value ? 1 : 0);
}

RecipeValue::RecipeValue(const char *value)
:node(std::make_shared<Node>())
{
  this->node->kind = Kind::String;
  this->node->strings.push_back(value);
}

RecipeValue::RecipeValue(const std::string &value)
:node(std::make_shared<Node>())
{
  this->node->kind = Kind::String;
  this->node->strings.push_back(value);
}

RecipeValue::RecipeValue(const std::vector<double> &values)
:node(std::make_shared<Node>())
{
  this->node->kind = Kind::Real;
  this->node->reals = values;
}

RecipeValue::RecipeValue(const std::vector<int> &values)
:node(std::make_shared<Node>())
{
  this->node->kind = Kind::Integer;
  this->node->logicals = values;
}

RecipeValue::RecipeValue(const std::vector<bool> &values)
:node(std::make_shared<Node>())
{
  this->node->kind = Kind::Logical;
  for (bool value : values)
    this->node->logicals.push_back(value ? 1 : 0);
}

RecipeValue::RecipeValue(const std::vector<std::string> &values)
:node(std::make_shared<Node>())
{
  this->node->kind = Kind::String;
  this->node->strings = values;
}

RecipeValue::RecipeValue(const arma::mat &value)
:node(std::make_shared<Node>())
{
  this->node->kind = Kind::Matrix;
  this->node->matrix = value;
}

RecipeValue::RecipeValue(const arma::colvec &value)
:node(std::make_shared<Node>())
{
  this->node->kind = Kind::Real;
  this->node->reals = arma::conv_to<std::vector<double>>::from(value);
}

RecipeValue::RecipeValue(const FunctionReference &value)
:node(std::make_shared<Node>())
{
  this->node->kind = Kind::Function;
  this->node->function = value;
}

RecipeValue::RecipeValue(const std::vector<std::vector<double>> &values)
:node(std::make_shared<Node>())
{
  this->node->kind = Kind::List;
  for (const std::vector<double> &value : values)
  {
    this->node->items.push_back(RecipeValue(value));
    this->node->keys.push_back("");
  }
}

void RecipeValue::set(const std::string &name, const RecipeValue &value)
{
  this->require(this->node->kind == Kind::List, "a list");
  Node &n = this->writable();
  for (size_t i = 0; i < n.keys.size(); ++i)
  {
    if (n.keys[i] == name)
    {
      n.items[i] = value;
      return;
    }
  }
  n.items.push_back(value);
  n.keys.push_back(name);
}

RecipeValue RecipeValue::list()
{
  auto n = std::make_shared<Node>();
  n->kind = Kind::List;
  return RecipeValue(n);
}

RecipeValue RecipeValue::list(std::initializer_list<RecipeValue> values)
{
  RecipeValue output = RecipeValue::list();
  for (const RecipeValue &value : values)
    output.push_back(value);
  return output;
}

RecipeValue RecipeValue::other(const std::string &description)
{
  auto n = std::make_shared<Node>();
  n->kind = Kind::Other;
  n->description = description;
  return RecipeValue(n);
}

RecipeValue RecipeValue::logicals(const std::vector<int> &values)
{
  auto n = std::make_shared<Node>();
  n->kind = Kind::Logical;
  n->logicals = values;
  return RecipeValue(n);
}

RecipeValue RecipeValue::integers(const std::vector<int> &values)
{
  auto n = std::make_shared<Node>();
  n->kind = Kind::Integer;
  n->logicals = values;
  return RecipeValue(n);
}

RecipeValue::Node &RecipeValue::writable()
{
  // copy on write: other values may share this node
  if (this->node.use_count() > 1)
    this->node = std::make_shared<Node>(*this->node);
  return *this->node;
}

void RecipeValue::push_back(const RecipeValue &value, const std::string &name)
{
  this->require(this->node->kind == Kind::List, "a list");
  Node &n = this->writable();
  n.items.push_back(value);
  n.keys.push_back(name);
}

RecipeValue::Kind RecipeValue::kind() const
{
  return this->node->kind;
}

std::string RecipeValue::kind_name() const
{
  switch (this->node->kind)
  {
    case Kind::Null: return "null";
    case Kind::Logical: return "logical";
    case Kind::Integer: return "integer";
    case Kind::Real: return "real";
    case Kind::String: return "string";
    case Kind::Matrix: return "matrix";
    case Kind::List: return "list";
    case Kind::Function: return "function";
    default: return "other (" + this->node->description + ")";
  }
}

void RecipeValue::require(bool condition, const std::string &wanted) const
{
  if (!condition)
    Rcpp::stop("Recipe: expected " + wanted + ", found " + this->kind_name() + ".");
}

bool RecipeValue::is_null() const { return this->node->kind == Kind::Null; }
bool RecipeValue::is_list() const { return this->node->kind == Kind::List; }
bool RecipeValue::is_function() const { return this->node->kind == Kind::Function; }
bool RecipeValue::is_string() const { return this->node->kind == Kind::String; }
bool RecipeValue::is_matrix() const { return this->node->kind == Kind::Matrix; }

bool RecipeValue::is_real() const
{
  return (this->node->kind == Kind::Real) || (this->node->kind == Kind::Matrix);
}

bool RecipeValue::is_numeric() const
{
  Kind k = this->node->kind;
  return (k == Kind::Real) || (k == Kind::Matrix) || (k == Kind::Integer) || (k == Kind::Logical);
}

bool RecipeValue::is_vector() const
{
  Kind k = this->node->kind;
  return (k == Kind::Real) || (k == Kind::Integer) || (k == Kind::Logical) ||
    (k == Kind::String) || (k == Kind::List) || (k == Kind::Matrix);
}

size_t RecipeValue::size() const
{
  switch (this->node->kind)
  {
    case Kind::Logical:
    case Kind::Integer: return this->node->logicals.size();
    case Kind::Real: return this->node->reals.size();
    case Kind::String: return this->node->strings.size();
    case Kind::Matrix: return this->node->matrix.n_elem;
    case Kind::List: return this->node->items.size();
    case Kind::Function: return 1;
    default: return 0;
  }
}

bool RecipeValue::contains(const std::string &name) const
{
  if (this->node->kind != Kind::List)
    return false;
  for (const std::string &key : this->node->keys)
  {
    if (key == name)
      return true;
  }
  return false;
}

const RecipeValue &RecipeValue::operator[](const std::string &name) const
{
  this->require(this->node->kind == Kind::List, "a list (looking up \"" + name + "\")");
  for (size_t i = 0; i < this->node->keys.size(); ++i)
  {
    if (this->node->keys[i] == name)
      return this->node->items[i];
  }
  Rcpp::stop("Recipe: no element named \"" + name + "\".");
}

const RecipeValue &RecipeValue::operator[](size_t index) const
{
  this->require(this->node->kind == Kind::List, "a list");
  if (index >= this->node->items.size())
    Rcpp::stop("Recipe: index " + std::to_string(index) + " out of range (list of length " +
               std::to_string(this->node->items.size()) + ").");
  return this->node->items[index];
}

const RecipeValue &RecipeValue::operator[](int index) const
{
  if (index < 0)
    Rcpp::stop("Recipe: negative index.");
  return (*this)[size_t(index)];
}

std::vector<std::string> RecipeValue::names() const
{
  if (this->node->kind != Kind::List)
    return std::vector<std::string>();
  return this->node->keys;
}

double RecipeValue::as_double() const
{
  switch (this->node->kind)
  {
    case Kind::Real:
      if (this->node->reals.empty()) break;
      return this->node->reals[0];
    case Kind::Integer:
    case Kind::Logical:
      if (this->node->logicals.empty()) break;
      if (this->node->logicals[0] == std::numeric_limits<int>::min())
        return std::numeric_limits<double>::quiet_NaN(); // R's NA
      return double(this->node->logicals[0]);
    case Kind::Matrix:
      if (this->node->matrix.n_elem == 0) break;
      return this->node->matrix(0);
    default:
      break;
  }
  Rcpp::stop("Recipe: expected a number, found " + this->kind_name() +
             (this->size() == 0 ? " of length 0" : "") + ".");
}

int RecipeValue::as_int() const
{
  return int(this->as_double());
}

size_t RecipeValue::as_size_t() const
{
  double value = this->as_double();
  if (value < 0.0)
    Rcpp::stop("Recipe: expected a non-negative whole number.");
  return size_t(value);
}

bool RecipeValue::as_bool() const
{
  return this->as_double() != 0.0;
}

std::string RecipeValue::as_string() const
{
  this->require((this->node->kind == Kind::String) && !this->node->strings.empty(), "a string");
  return this->node->strings[0];
}

arma::mat RecipeValue::as_matrix() const
{
  if (this->node->kind == Kind::Matrix)
    return this->node->matrix;
  std::vector<double> values = this->as_doubles();
  return arma::mat(values.data(), values.size(), values.empty() ? 0 : 1);
}

arma::colvec RecipeValue::as_colvec() const
{
  return arma::vectorise(this->as_matrix());
}

arma::rowvec RecipeValue::as_rowvec() const
{
  return arma::vectorise(this->as_matrix()).t();
}

std::vector<double> RecipeValue::as_doubles() const
{
  switch (this->node->kind)
  {
    case Kind::Real:
      return this->node->reals;
    case Kind::Integer:
    case Kind::Logical:
    {
      std::vector<double> output;
      for (int value : this->node->logicals)
        output.push_back(double(value));
      return output;
    }
    case Kind::Matrix:
      return arma::conv_to<std::vector<double>>::from(arma::vectorise(this->node->matrix));
    case Kind::Null:
      // as R's as.numeric(NULL) and Rcpp's as<std::vector<double>>(R_NilValue)
      // (e.g. order_of_mcmc in a model without MCMC moves, read by EnKI)
      return std::vector<double>();
    default:
      break;
  }
  Rcpp::stop("Recipe: expected numbers, found " + this->kind_name() + ".");
}

std::vector<int> RecipeValue::as_ints() const
{
  std::vector<int> output;
  for (double value : this->as_doubles())
    output.push_back(int(value));
  return output;
}

std::vector<size_t> RecipeValue::as_size_ts() const
{
  std::vector<size_t> output;
  for (double value : this->as_doubles())
  {
    if (value < 0.0)
      Rcpp::stop("Recipe: expected non-negative whole numbers.");
    output.push_back(size_t(value));
  }
  return output;
}

std::vector<std::string> RecipeValue::as_strings() const
{
  this->require(this->node->kind == Kind::String, "strings");
  return this->node->strings;
}

std::vector<std::vector<double>> RecipeValue::as_list_of_doubles() const
{
  this->require(this->node->kind == Kind::List, "a list of numeric vectors");
  std::vector<std::vector<double>> output;
  for (const RecipeValue &item : this->node->items)
    output.push_back(item.as_doubles());
  return output;
}

const FunctionReference &RecipeValue::as_function() const
{
  this->require(this->node->kind == Kind::Function, "a function");
  return this->node->function;
}

static void json_string(std::string &out, const std::string &value)
{
  out += '"';
  for (char c : value)
  {
    switch (c)
    {
      case '"': out += "\\\""; break;
      case '\\': out += "\\\\"; break;
      case '\n': out += "\\n"; break;
      case '\t': out += "\\t"; break;
      case '\r': out += "\\r"; break;
      default:
        if (static_cast<unsigned char>(c) < 0x20)
        {
          char buffer[8];
          std::snprintf(buffer, sizeof(buffer), "\\u%04x", c);
          out += buffer;
        }
        else
        {
          out += c;
        }
    }
  }
  out += '"';
}

static void json_number(std::string &out, double value)
{
  if (std::isnan(value)) { out += "\"NaN\""; return; }
  if (std::isinf(value)) { out += (value > 0 ? "\"Inf\"" : "\"-Inf\""); return; }
  char buffer[32];
  std::snprintf(buffer, sizeof(buffer), "%.17g", value);
  out += buffer;
}

void RecipeValue::write_json(std::string &out) const
{
  const Node &n = *this->node;
  switch (n.kind)
  {
    case Kind::Null:
      out += "null";
      return;
    case Kind::Logical:
    case Kind::Integer:
    {
      out += (n.kind == Kind::Logical ? "{\"logical\":[" : "{\"integer\":[");
      for (size_t i = 0; i < n.logicals.size(); ++i)
      {
        if (i > 0) out += ',';
        if (n.logicals[i] == std::numeric_limits<int>::min())
          out += "\"NA\"";
        else if (n.kind == Kind::Logical)
          out += (n.logicals[i] ? "true" : "false");
        else
          out += std::to_string(n.logicals[i]);
      }
      out += "]}";
      return;
    }
    case Kind::Real:
    {
      out += "{\"real\":[";
      for (size_t i = 0; i < n.reals.size(); ++i)
      {
        if (i > 0) out += ',';
        json_number(out, n.reals[i]);
      }
      out += "]}";
      return;
    }
    case Kind::String:
    {
      out += "{\"string\":[";
      for (size_t i = 0; i < n.strings.size(); ++i)
      {
        if (i > 0) out += ',';
        json_string(out, n.strings[i]);
      }
      out += "]}";
      return;
    }
    case Kind::Matrix:
    {
      out += "{\"matrix\":{\"rows\":" + std::to_string(n.matrix.n_rows) +
        ",\"cols\":" + std::to_string(n.matrix.n_cols) + ",\"values\":[";
      for (size_t i = 0; i < n.matrix.n_elem; ++i)
      {
        if (i > 0) out += ',';
        json_number(out, n.matrix(i));
      }
      out += "]}}";
      return;
    }
    case Kind::List:
    {
      out += "{\"list\":[";
      for (size_t i = 0; i < n.items.size(); ++i)
      {
        if (i > 0) out += ',';
        out += "{\"name\":";
        json_string(out, n.keys[i]);
        out += ",\"value\":";
        n.items[i].write_json(out);
        out += '}';
      }
      out += "]}";
      return;
    }
    case Kind::Function:
    {
      out += "{\"function\":{\"language\":";
      json_string(out, n.function.language);
      out += ",\"name\":";
      json_string(out, n.function.name);
      out += "}}";
      return;
    }
    default:
      out += "{\"other\":";
      json_string(out, n.description);
      out += '}';
  }
}

std::string RecipeValue::to_json() const
{
  std::string out;
  this->write_json(out);
  return out;
}

}
