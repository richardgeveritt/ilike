// Lowering the ilike language (dev/ilike-language.md; dev/plan.md 3.14, I3b):
// from the syntax trees of the files of a recipe to the recipe the engine
// reads — the same tree the old parser makes — with each user function
// marked by a placeholder {ilike_function, role} for the binding to compile
// and replace by the compiled function (role: the old block name, which fixes
// the function's signature). Small functions that the old format had users
// write (the name of a data variable, a list of factor numbers, literal data)
// are generated here as C++ code.

#include <RcppArmadillo.h>

#include <cctype>
#include <cmath>
#include <map>
#include <set>
#include <sstream>

#include "ilike_language.h"
#include "builtin_registry.h"

namespace ilike
{

namespace
{

std::string number_text(double value)
{
  if (std::isinf(value))
    return value > 0 ? "std::numeric_limits<double>::infinity()" : "-std::numeric_limits<double>::infinity()";
  std::ostringstream out;
  out.precision(17);
  out << value;
  return out.str();
}

std::string quoted(const std::string &text)
{
  std::string output = "\"";
  for (char c : text)
  {
    if (c == '"' || c == '\\')
      output.push_back('\\');
    output.push_back(c);
  }
  return output + "\"";
}

bool role_exists(const std::string &name, BuiltinRole role)
{
  for (const Builtin &entry : builtin_registry())
  {
    std::string a = entry.name, b = name;
    if (a.size() != b.size())
      continue;
    bool same = true;
    for (size_t i = 0; i < a.size(); ++i)
      if (std::tolower((unsigned char)a[i]) != std::tolower((unsigned char)b[i]))
        same = false;
    if (same && entry.has_role(role))
      return true;
  }
  return false;
}

class Lowering
{

public:

  explicit Lowering(const RecipeValue &trees_in) : trees(trees_in)
  {
    this->recipe = RecipeValue::list();
    this->data = RecipeValue::list();
    this->algorithms = RecipeValue::list();
    this->constants = RecipeValue::list();
  }

  RecipeValue run()
  {
    // the C++ code first, so that references to it can be checked
    for (size_t t = 0; t < this->trees.size(); ++t)
    {
      const RecipeValue &sections = this->trees[t]["sections"];
      for (size_t s = 0; s < sections.size(); ++s)
      {
        const RecipeValue &section = sections[s];
        std::string kind = section["section"].as_string();
        if (kind == "cpp")
          this->user_code += section["code"].as_string() + "\n";
      }
    }

    // the imports, so that the model may refer to them wherever they are
    for (size_t t = 0; t < this->trees.size(); ++t)
    {
      this->file = this->trees[t]["file"].as_string();
      const RecipeValue &sections = this->trees[t]["sections"];
      for (size_t s = 0; s < sections.size(); ++s)
        if (sections[s]["section"].as_string() == "import")
        {
          std::string alias = sections[s]["alias"].as_string();
          std::string path = sections[s]["path"].as_string();
          if (this->modules.count(alias) || (this->foreign_languages.count(alias) && alias != "r" && alias != "py" && alias != "jl"))
            this->error(sections[s], "\"" + alias + "\" is imported more than once.");
          std::string extension = path.substr(path.find_last_of('.') == std::string::npos ? path.size() : path.find_last_of('.') + 1);
          std::string language = (extension == "R" || extension == "r") ? "r" : (extension == "py") ? "py" : (extension == "jl") ? "jl" : "";
          if (!language.empty())
          {
            // a file of functions in another language, used as alias::f
            this->foreign_languages[alias] = language;
            RecipeValue file_entry = RecipeValue::list();
            file_entry.push_back(RecipeValue(language), "language");
            file_entry.push_back(RecipeValue(path), "path");
            file_entry.push_back(RecipeValue(this->file), "importer");
            this->foreign_files.push_back(file_entry);
          }
          else
            this->modules[alias] = std::make_pair(path, this->file);
        }
    }

    // the model (factors and their names) before the algorithms, which
    // refer to factors by name
    // (constants and data first: calls to other languages need to know which
    // names are data and constants; then the model; then the algorithms)
    for (int pass = 0; pass < 3; ++pass)
    {
      for (size_t t = 0; t < this->trees.size(); ++t)
      {
        this->file = this->trees[t]["file"].as_string();
        const RecipeValue &sections = this->trees[t]["sections"];
        for (size_t s = 0; s < sections.size(); ++s)
        {
          const RecipeValue &section = sections[s];
          std::string kind = section["section"].as_string();
          if (pass == 0)
          {
            if (kind == "data")
              this->data_section(section);
            else if (kind == "constants")
              this->constants_section(section);
          }
          else if (pass == 1)
          {
            if (kind == "model")
            {
              if (!this->generation_decided)
                this->decide_generation();
              this->model_section(section);
            }
            else if (kind == "r" || kind == "py" || kind == "jl")
            {
              RecipeValue code_entry = RecipeValue::list();
              code_entry.push_back(RecipeValue(kind), "language");
              code_entry.push_back(section["code"], "code");
              this->foreign_code.push_back(code_entry);
            }
            else if (kind == "import")
              continue; // (read first)
          }
          else if (kind == "algorithm")
            this->algorithm_section(section);
        }
      }
      if (pass == 1)
        this->finish_model_statements();
    }
    this->finish_data();

    if (this->data.size() > 0)
      this->recipe.push_back(this->data, "data");
    if (!this->factor_entries.empty())
    {
      RecipeValue factors = RecipeValue::list();
      for (const RecipeValue &factor : this->factor_entries)
        factors.push_back(factor);
      this->recipe.push_back(factors, "factor");
    }
    if (this->transitions.size() > 0)
      this->recipe.push_back(this->transitions, "transition_model");
    if (this->model_methods.size() > 0)
      this->recipe.push_back(this->model_methods, "method");

    RecipeValue output = RecipeValue::list();
    output.push_back(this->recipe, "recipe");
    output.push_back(this->algorithms, "algorithms");
    output.push_back(this->constants, "constants");
    output.push_back(RecipeValue(this->declared_constants), "declared_constants");
    output.push_back(RecipeValue(this->user_code), "user_code");
    std::string all_generated = this->generated_code;
    if (!this->model_code.empty())
      all_generated += this->model_header + this->model_code + this->deferred_node_code + this->deferred_function_code;
    output.push_back(RecipeValue(all_generated), "generated_code");
    output.push_back(RecipeValue(std::vector<std::string>(this->model_constants_used.begin(), this->model_constants_used.end())),
                     "model_constants");
    output.push_back(this->foreign_code, "foreign_code");
    output.push_back(this->foreign_files, "foreign_files");
    output.push_back(RecipeValue(this->uses_foreign), "uses_foreign");
    return output;
  }

private:

  RecipeValue trees;
  std::string file;
  RecipeValue recipe;
  std::vector<RecipeValue> factor_entries;
  RecipeValue data;
  RecipeValue algorithms;
  RecipeValue constants;
  std::vector<std::string> declared_constants;
  std::string user_code;
  std::string generated_code;
  std::map<std::string, size_t> factor_numbers;   // name -> 1, 2, ...
  std::set<std::string> defined_names;
  std::vector<std::pair<std::string, RecipeValue>> literal_data;
  RecipeValue transitions = RecipeValue::list();   // transition_model
  RecipeValue model_methods = RecipeValue::list(); // method entries of the model (ssm)
  std::map<std::string, std::pair<std::string, std::string>> modules; // alias -> (path, importing file)
  // functions in other languages (I5): namespace -> language, the code
  // sections and imported files to load, the names known to be data or
  // constants (for the translation of calls into R expressions)
  std::map<std::string, std::string> foreign_languages = {{"r", "r"}, {"py", "py"}, {"jl", "jl"}};
  RecipeValue foreign_code = RecipeValue::list();
  RecipeValue foreign_files = RecipeValue::list();
  std::set<std::string> data_names;
  std::set<std::string> constant_names;
  bool uses_foreign = false;

  // model statements (I4a)
  struct ModelLoop { std::string index; RecipeValue from; RecipeValue to; };
  struct ModelStatement
  {
    std::string op;            // "~" or "<-"
    std::string variable;      // the target's name
    RecipeValue indices;       // the target's index expressions (a list; empty for none)
    RecipeValue value;         // the distribution call, or the expression
    std::vector<ModelLoop> loops;
    RecipeValue node;          // for messages
    bool observed = false;     // the target is data
  };
  std::vector<ModelStatement> model_statements;
  bool generate = false;
  bool generation_decided = false;
  std::map<std::string, size_t> model_factor_position; // "model_prior" or a data variable -> factor entry
  std::set<std::string> model_constants_used;
  std::string model_header;            // static literal data
  std::string model_code;              // helpers and declarations
  std::string deferred_node_code;      // deterministic nodes
  std::string deferred_function_code;  // the factors' functions
  bool data_function = false;
  size_t generated_count = 0;

  // --- messages ------------------------------------------------------------------

  [[noreturn]] void error(const RecipeValue &node, const std::string &message) const
  {
    std::string where = this->file;
    if (node.is_list() && node.contains("line"))
      where = where + ":" + std::to_string(node["line"].as_int()) + ":" + std::to_string(node["column"].as_int());
    Rcpp::stop(where + ": " + message);
  }

  // --- expressions ---------------------------------------------------------------

  // loops are model statements: the parser reads them anywhere, the other sections refuse them here
  void refuse_loop(const RecipeValue &statement, const std::string &section) const
  {
    if (statement.contains("keyword") && statement["keyword"].as_string() == "for")
      this->error(statement, "loops (for) are only allowed in the model section, not in " + section + ".");
  }

  static std::string kind_of(const RecipeValue &expression)
  {
    return expression.is_list() && expression.contains("expression") ? expression["expression"].as_string() : "";
  }

  bool is_function_reference(const RecipeValue &expression) const
  {
    return kind_of(expression) == "name" && !expression["namespace"].as_string().empty() &&
           expression["namespace"].as_string() != "ilike";
  }

  // the name of a C++ function referred to as cpp::f (checked to be defined)
  std::string cpp_function(const RecipeValue &expression, const std::string &what)
  {
    if (kind_of(expression) != "name" || expression["namespace"].as_string().empty())
      this->error(expression, what + " must be a function, written cpp::name.");
    std::string name_space = expression["namespace"].as_string();
    if (name_space != "cpp")
      this->error(expression, "unknown namespace \"" + name_space + "\" (functions are written cpp::name).");
    std::string name = expression["name"].as_string();
    if (find_cpp_function_signature(this->user_code, name).empty())
      this->error(expression, "cpp::" + name + " is not defined in a cpp section.");
    return name;
  }

  RecipeValue placeholder(const std::string &function, const std::string &role) const
  {
    RecipeValue output = RecipeValue::list();
    output.push_back(RecipeValue(function), "ilike_function");
    output.push_back(RecipeValue(role), "role");
    return output;
  }

  // --- functions in other languages (I5) ------------------------------------------

  // TRUE for r::f, py::f, jl::f, alias::f (of an imported file), or a call to one
  bool is_foreign(const RecipeValue &expression) const
  {
    std::string kind = kind_of(expression);
    if (kind != "name" && kind != "call")
      return false;
    std::string name_space = expression["namespace"].as_string();
    return !name_space.empty() && name_space != "cpp" && name_space != "ilike" && this->foreign_languages.count(name_space);
  }

  // a function: cpp::f, or a function (or call) in another language
  bool is_function_value(const RecipeValue &expression) const
  {
    return this->is_function_reference(expression) || this->is_foreign(expression);
  }

  // the arguments with which a function in another language is called when
  // it is given without them (r::f): the simplest form of each role, as in
  // the old format's R expressions
  static std::string default_arguments(const std::string &role)
  {
    static const std::map<std::string, std::string> table = {
      {"data", ""}, {"evaluate_log_prior", "parameters"}, {"simulate_prior", ""},
      {"evaluate_gradient_log_prior", "variable, parameters"},
      {"evaluate_second_gradient_log_prior", "variable1, variable2, parameters"},
      {"evaluate_log_likelihood", "parameters, data"},
      {"evaluate_gradient_log_likelihood", "variable, parameters, data"},
      {"evaluate_second_gradient_log_likelihood", "variable1, variable2, parameters, data"},
      {"simulate_data_model", "parameters"}, {"summary_statistics", "data"},
      {"simulate_m_proposal", "parameters"}, {"simulate_mh_proposal", "parameters"},
      {"simulate_unadjusted_proposal", "parameters"},
      {"evaluate_log_mh_proposal", "proposed_parameters, parameters"},
      {"simulate_imh_proposal", ""}, {"evaluate_log_imh_proposal", "proposed_parameters"},
      {"simulate_is_proposal", ""}, {"evaluate_log_is_proposal", "proposed_parameters"},
      {"simulate_transition_model", "parameters"},
      {"evaluate_log_transition_model", "proposed_parameters, parameters"},
      {"simulate_transition_proposal", "parameters, data"},
      {"evaluate_log_transition_proposal", "proposed_parameters, parameters, data"},
      {"nonlinear_gaussian_data_function", "parameters"}, {"nonlinear_gaussian_transition_function", "parameters"},
      {"linear_gaussian_data_matrix", "parameters"}, {"linear_gaussian_data_covariance", "parameters"},
      {"linear_gaussian_transition_matrix", "parameters"}, {"linear_gaussian_transition_covariance", "parameters"},
      {"nonlinear_gaussian_data_covariance", "parameters"}, {"nonlinear_gaussian_transition_covariance", "parameters"},
      {"mcmc_weights", ""}};
    auto found = table.find(role);
    return found == table.end() ? std::string("") : found->second;
  }

  static std::string r_number(double value)
  {
    if (std::isnan(value)) return "NaN";
    if (std::isinf(value)) return value > 0 ? "Inf" : "-Inf";
    std::ostringstream out;
    out.precision(17);
    out << value;
    return out.str();
  }

  // an expression as R code, for the R binding (which runs Python and Julia
  // functions through R): a bare name is a parameter (parameters$name),
  // unless it is a data variable (data$name), a constant, or one of the
  // objects a role passes (parameters, data, proposed_parameters,
  // proposal_parameters, variable, variable1, variable2); calls inside are
  // calls in that language
  std::string r_text(const RecipeValue &expression) const
  {
    static const std::set<std::string> special = {"parameters", "data", "proposed_parameters", "proposal_parameters",
                                                  "variable", "variable1", "variable2"};
    std::string kind = kind_of(expression);
    if (kind == "number")
    {
      if (expression.contains("integer"))
        return r_number(expression["value"].as_double()) + "L";
      return r_number(expression["value"].as_double());
    }
    if (kind == "string")
      return quoted(expression["value"].as_string());
    if (kind == "logical")
      return expression["value"].as_bool() ? "TRUE" : "FALSE";
    if (kind == "name")
    {
      std::string name = expression["name"].as_string();
      if (!expression["namespace"].as_string().empty())
        return name;
      if (special.count(name) || this->constant_names.count(name))
        return name;
      if (this->data_names.count(name))
        return "data$" + name;
      return "parameters$" + name;
    }
    if (kind == "member")
      return this->r_text(expression["target"]) + "$" + expression["name"].as_string();
    if (kind == "call")
    {
      const RecipeValue &arguments = expression["arguments"];
      std::vector<std::string> keys = arguments.names();
      std::string text = expression["function"].as_string() + "(";
      for (size_t k = 0; k < arguments.size(); ++k)
        text += (k ? ", " : "") + ((k < keys.size() && !keys[k].empty()) ? keys[k] + " = " : std::string("")) +
                this->r_text(arguments[k]);
      return text + ")";
    }
    if (kind == "vector")
    {
      const RecipeValue &elements = expression["elements"];
      std::string text = "c(";
      for (size_t k = 0; k < elements.size(); ++k)
        text += (k ? ", " : "") + this->r_text(elements[k]);
      return text + ")";
    }
    if (kind == "index")
    {
      const RecipeValue &indices = expression["indices"];
      std::string text = this->r_text(expression["target"]) + "[";
      for (size_t k = 0; k < indices.size(); ++k)
        text += (k ? ", " : "") + this->r_text(indices[k]);
      return text + "]";
    }
    if (kind == "binary")
      return "(" + this->r_text(expression["left"]) + " " + expression["operator"].as_string() + " " +
             this->r_text(expression["right"]) + ")";
    if (kind == "unary")
      return "(-" + this->r_text(expression["operand"]) + ")";
    this->error(expression, "this cannot be passed to a function in another language.");
  }

  // a function in another language used in a role: {ilike_foreign: language,
  // expression (R code calling it), role, output (the variable a data
  // function gives, or "")}; the R binding generates its C++ wrapper
  RecipeValue foreign_placeholder(const RecipeValue &expression, const std::string &role, const std::string &output)
  {
    if (role.compare(0, 6, "adapt_") == 0 || role.find("factor_index") != std::string::npos)
      this->error(expression, "this function must be written in C++ (cpp::f).");
    std::string language = this->foreign_languages.at(expression["namespace"].as_string());
    std::string text;
    if (kind_of(expression) == "name")
      text = expression["name"].as_string() + "(" + default_arguments(role) + ")";
    else
      text = this->r_text(expression);
    this->uses_foreign = true;
    RecipeValue output_value = RecipeValue::list();
    output_value.push_back(RecipeValue(language), "ilike_foreign");
    output_value.push_back(RecipeValue(text), "expression");
    output_value.push_back(RecipeValue(role), "role");
    output_value.push_back(RecipeValue(output), "output");
    return output_value;
  }

  // the placeholder of a function used in a role: C++ or another language
  RecipeValue function_placeholder(const RecipeValue &expression, const std::string &role, const std::string &what)
  {
    if (this->is_foreign(expression))
      return this->foreign_placeholder(expression, role, "");
    return this->placeholder(this->cpp_function(expression, what), role);
  }

  // the name of a built-in call (ilike:: accepted as a synonym)
  std::string builtin_name(const RecipeValue &call) const
  {
    std::string name_space = call["namespace"].as_string();
    if (!name_space.empty() && name_space != "ilike")
      this->error(call, name_space + "::" + call["function"].as_string() + " is not a built-in (built-ins are written without a prefix, or with ilike::).");
    return call["function"].as_string();
  }

  // a constant numeric expression, folded; false if it is not one
  bool numeric(const RecipeValue &expression, arma::mat &value) const
  {
    std::string kind = kind_of(expression);
    if (kind == "number")
    {
      value = arma::mat(1, 1, arma::fill::value(expression["value"].as_double()));
      return true;
    }
    if (kind == "unary")
    {
      if (!this->numeric(expression["operand"], value))
        return false;
      value = -value;
      return true;
    }
    if (kind == "binary")
    {
      std::string op = expression["operator"].as_string();
      arma::mat left, right;
      if (!this->numeric(expression["left"], left) || !this->numeric(expression["right"], right))
        return false;
      if (op == ":")
      {
        if (left.n_elem != 1 || right.n_elem != 1)
          return false;
        double a = left(0), b = right(0);
        std::vector<double> values;
        if (a <= b)
          for (double x = a; x <= b + 1e-10; x += 1.0) values.push_back(x);
        else
          for (double x = a; x >= b - 1e-10; x -= 1.0) values.push_back(x);
        value = arma::mat(values.data(), values.size(), 1);
        return true;
      }
      bool scalar_left = (left.n_elem == 1), scalar_right = (right.n_elem == 1);
      if (!scalar_left && !scalar_right && arma::size(left) != arma::size(right))
        this->error(expression, "the sizes of the two sides of \"" + op + "\" differ.");
      auto apply = [&](double a, double b) {
        if (op == "+") return a + b;
        if (op == "-") return a - b;
        if (op == "*") return a * b;
        if (op == "/") return a / b;
        return std::pow(a, b);
      };
      arma::mat shape = scalar_left ? right : left;
      value = arma::mat(arma::size(shape));
      for (arma::uword i = 0; i < value.n_elem; ++i)
        value(i) = apply(scalar_left ? left(0) : left(i), scalar_right ? right(0) : right(i));
      return true;
    }
    if (kind == "vector")
    {
      const RecipeValue &elements = expression["elements"];
      std::vector<arma::mat> items;
      for (size_t i = 0; i < elements.size(); ++i)
      {
        arma::mat item;
        if (!this->numeric(elements[i], item))
          return false;
        items.push_back(item);
      }
      if (items.empty())
      {
        value = arma::mat(0, 1);
        return true;
      }
      bool all_scalars = true, all_rows = true;
      for (const arma::mat &item : items)
      {
        if (item.n_elem != 1) all_scalars = false;
        if (!(item.n_cols == 1 && item.n_rows == items[0].n_rows && kind_of(elements[0]) == "vector")) all_rows = false;
      }
      if (all_scalars)
      {
        value = arma::mat(items.size(), 1);
        for (size_t i = 0; i < items.size(); ++i) value(i) = items[i](0);
        return true;
      }
      if (all_rows)
      {
        // [[a, b], [c, d]]: a matrix given by rows
        value = arma::mat(items.size(), items[0].n_rows);
        for (size_t i = 0; i < items.size(); ++i) value.row(i) = items[i].t();
        return true;
      }
      this->error(expression, "a vector of numbers, or a matrix given as rows of equal length, is needed.");
    }
    if (kind == "call" && call_is_helper(expression))
      return this->helper(expression, value);
    return false;
  }

  static bool call_is_helper(const RecipeValue &call)
  {
    std::string f = call["function"].as_string();
    return call["namespace"].as_string().empty() && (f == "seq" || f == "rep" || f == "diag");
  }

  // seq(from, to, length = n | by = b), rep(x, n), diag(n) or diag(v)
  bool helper(const RecipeValue &call, arma::mat &value) const
  {
    std::string f = call["function"].as_string();
    const RecipeValue &arguments = call["arguments"];
    std::vector<std::string> keys = arguments.names();
    std::map<std::string, arma::mat> named;
    std::vector<arma::mat> positional;
    for (size_t i = 0; i < arguments.size(); ++i)
    {
      arma::mat v;
      if (!this->numeric(arguments[i], v))
        return false;
      if (i < keys.size() && !keys[i].empty()) named[keys[i]] = v; else positional.push_back(v);
    }
    if (f == "seq")
    {
      if (positional.size() < 2)
        this->error(call, "seq(from, to, length = n) or seq(from, to, by = b).");
      double from = positional[0](0), to = positional[1](0);
      std::vector<double> values;
      if (named.count("length") || named.count("length.out"))
      {
        double n = named.count("length") ? named["length"](0) : named["length.out"](0);
        if (n < 1 || n != std::floor(n))
          this->error(call, "seq: length must be a positive whole number.");
        for (double k = 0; k < n; ++k)
          values.push_back(n == 1 ? from : from + (to - from) * k / (n - 1));
      }
      else
      {
        double by = named.count("by") ? named["by"](0) : (positional.size() > 2 ? positional[2](0) : (to >= from ? 1.0 : -1.0));
        if (by == 0 || (to - from) * by < 0)
          this->error(call, "seq: the step does not lead from \"from\" to \"to\".");
        for (double x = from; (by > 0) ? x <= to + 1e-10 * std::fabs(by) : x >= to - 1e-10 * std::fabs(by); x += by)
          values.push_back(x);
      }
      value = arma::mat(values.data(), values.size(), 1);
      return true;
    }
    if (f == "rep")
    {
      if (positional.size() != 2 || positional[1].n_elem != 1)
        this->error(call, "rep(x, n).");
      value = arma::repmat(arma::vectorise(positional[0]), size_t(positional[1](0)), 1);
      return true;
    }
    // diag
    if (positional.size() != 1)
      this->error(call, "diag(n) or diag(v).");
    if (positional[0].n_elem == 1)
      value = arma::eye(size_t(positional[0](0)), size_t(positional[0](0)));
    else
      value = arma::diagmat(arma::vectorise(positional[0]));
    return true;
  }

  // the value of an argument: a number, vector or matrix (folded), a
  // logical, a string, or a name (of a constant or a variable)
  RecipeValue value_of(const RecipeValue &expression)
  {
    std::string kind = kind_of(expression);
    if (kind == "string")
      return RecipeValue(expression["value"].as_string());
    if (kind == "logical")
      return RecipeValue::logicals(std::vector<int>(1, expression["value"].as_bool() ? 1 : 0));
    if (kind == "name")
    {
      if (this->is_function_reference(expression))
        this->error(expression, "a function (" + expression["namespace"].as_string() + "::" +
                    expression["name"].as_string() + ") cannot be used here.");
      return RecipeValue(expression["name"].as_string());
    }
    arma::mat value;
    if (this->numeric(expression, value))
    {
      if (value.n_elem == 1 && kind != "vector")
        return RecipeValue(value(0));
      if (value.n_cols == 1 && !(kind == "vector" && kind_of(expression["elements"][0]) == "vector") &&
          !(kind == "call" && expression["function"].as_string() == "diag"))
        return RecipeValue(arma::conv_to<std::vector<double>>::from(arma::vectorise(value)));
      return RecipeValue(value);
    }
    if (kind == "binary" || kind == "unary")
      this->error(expression, "arithmetic is only available on numbers for now (not on constants or variables).");
    if (kind == "call")
      this->error(expression, "a call to " + expression["function"].as_string() + " cannot be used as a value here.");
    this->error(expression, "a value is needed here.");
  }

  // a built-in call as the old recipe's node {type, variables, parameters}:
  // the variable is the given one, or the first positional argument (or the
  // argument "variable"); with_variable = false for settings without one
  RecipeValue builtin_node(const RecipeValue &call, bool with_variable, const std::string &variable = "")
  {
    if (kind_of(call) != "call")
      this->error(call, "a built-in call is needed here (e.g. norm_rw(theta, 1)).");
    std::string type = this->builtin_name(call);
    const RecipeValue &arguments = call["arguments"];
    std::vector<std::string> keys = arguments.names();
    std::string the_variable = variable;
    bool variable_taken = !variable.empty() || !with_variable;
    RecipeValue parameters = RecipeValue::list();
    for (size_t i = 0; i < arguments.size(); ++i)
    {
      std::string key = (i < keys.size()) ? keys[i] : "";
      if (with_variable && variable.empty() && (key == "variable" || (key.empty() && !variable_taken)))
      {
        if (kind_of(arguments[i]) != "name" || this->is_function_reference(arguments[i]))
          this->error(arguments[i], type + ": the first argument is the name of a variable.");
        the_variable = arguments[i]["name"].as_string();
        variable_taken = true;
        continue;
      }
      parameters.push_back(this->value_of(arguments[i]), key);
    }
    if (with_variable && the_variable.empty())
      this->error(call, type + ": the variable is missing.");
    RecipeValue node = RecipeValue::list();
    node.push_back(RecipeValue(type), "type");
    if (with_variable)
      node.push_back(RecipeValue(the_variable), "variables");
    node.push_back(parameters, "parameters");
    return node;
  }

  // keyword arguments of a call by name, positional ones in order
  struct Arguments
  {
    std::vector<RecipeValue> positional;
    std::map<std::string, RecipeValue> named;
    bool has(const std::string &name) const { return named.count(name) > 0; }
  };

  Arguments split_arguments(const RecipeValue &call) const
  {
    Arguments output;
    const RecipeValue &arguments = call["arguments"];
    std::vector<std::string> keys = arguments.names();
    for (size_t i = 0; i < arguments.size(); ++i)
    {
      std::string key = (i < keys.size()) ? keys[i] : "";
      if (key.empty())
        output.positional.push_back(arguments[i]);
      else
        output.named[key] = arguments[i];
    }
    return output;
  }

  void allow_only(const RecipeValue &call, const Arguments &arguments, const std::set<std::string> &allowed,
                  size_t positional_allowed) const
  {
    for (const auto &kv : arguments.named)
      if (!allowed.count(kv.first))
      {
        std::string list;
        for (const std::string &a : allowed) list = list + (list.empty() ? "" : ", ") + a;
        this->error(kv.second, call["function"].as_string() + ": no argument called " + kv.first +
                    " (it takes " + list + ").");
      }
    if (arguments.positional.size() > positional_allowed)
      this->error(call, call["function"].as_string() + ": too many arguments without names.");
  }

  // --- generated functions ---------------------------------------------------------

  std::string generate_string_function(const std::string &value)
  {
    std::string name = "ilike_generated_string_" + std::to_string(++this->generated_count);
    this->generated_code += "std::string " + name + "()\n{\n  return " + quoted(value) + ";\n}\n\n";
    return name;
  }

  // [names] (or numbers) of factors as their numbers 1, 2, ...
  RecipeValue factor_numbers_of(const RecipeValue &expression)
  {
    arma::mat numbers;
    if (this->numeric(expression, numbers))
      return RecipeValue(arma::conv_to<std::vector<double>>::from(arma::vectorise(numbers)));
    std::vector<double> output;
    std::vector<RecipeValue> items;
    if (kind_of(expression) == "vector")
      for (size_t i = 0; i < expression["elements"].size(); ++i) items.push_back(expression["elements"][i]);
    else
      items.push_back(expression);
    for (const RecipeValue &item : items)
    {
      if (kind_of(item) != "name")
        this->error(item, "a list of factor names is needed (e.g. [lik, theta]).");
      auto found = this->factor_numbers.find(item["name"].as_string());
      if (found == this->factor_numbers.end())
        this->error(item, "no factor called " + item["name"].as_string() + ".");
      output.push_back(double(found->second));
    }
    return RecipeValue(output);
  }

  std::string generate_factor_index(const RecipeValue &names_expression)
  {
    std::vector<std::string> names;
    if (kind_of(names_expression) == "name")
      names.push_back(names_expression["name"].as_string());
    else if (kind_of(names_expression) == "vector")
    {
      const RecipeValue &elements = names_expression["elements"];
      for (size_t i = 0; i < elements.size(); ++i)
      {
        if (kind_of(elements[i]) != "name")
          this->error(elements[i], "a list of factor names is needed (e.g. [lik, theta]).");
        names.push_back(elements[i]["name"].as_string());
      }
    }
    else
      this->error(names_expression, "a list of factor names is needed (e.g. [lik, theta]).");
    std::string numbers;
    for (const std::string &n : names)
    {
      auto found = this->factor_numbers.find(n);
      if (found == this->factor_numbers.end())
      {
        std::string list;
        for (const auto &kv : this->factor_numbers) list = list + (list.empty() ? "" : ", ") + kv.first;
        this->error(names_expression, "no factor called " + n + " (factors: " + (list.empty() ? "none" : list) + ").");
      }
      numbers = numbers + (numbers.empty() ? "" : ", ") + std::to_string(found->second);
    }
    std::string name = "ilike_generated_factor_index_" + std::to_string(++this->generated_count);
    this->generated_code += "NumericVector " + name + "()\n{\n  return NumericVector::create(" + numbers + ");\n}\n\n";
    return name;
  }

  // --- sections --------------------------------------------------------------------

  void define(const RecipeValue &statement, const std::string &name)
  {
    if (this->defined_names.count(name))
      this->error(statement, "\"" + name + "\" is defined more than once.");
    this->defined_names.insert(name);
  }

  void add_factor(const RecipeValue &statement, const std::string &name, const RecipeValue &factor)
  {
    this->define(statement, name);
    this->factor_entries.push_back(factor);
    this->factor_numbers[name] = this->factor_entries.size();
  }

  void constants_section(const RecipeValue &section)
  {
    // ("constant": a constant of unstated type, as convert_ilike_file() writes
    // the old format's p1, p2, ...)
    static const std::set<std::string> types = {"real", "int", "integer", "logical", "string", "vector", "matrix", "constant"};
    const RecipeValue &statements = section["statements"];
    for (size_t i = 0; i < statements.size(); ++i)
    {
      const RecipeValue &statement = statements[i];
      this->refuse_loop(statement, "constants");
      std::string keyword = statement["keyword"].as_string();
      if (keyword.empty() && statement["operator"].as_string() != "=")
        this->error(statement, "a constant is declared as \"type name\" or \"type name = default\" (types: real, int, logical, string, vector, matrix, constant).");
      if (!keyword.empty() && !types.count(keyword))
        this->error(statement, "unknown type \"" + keyword + "\" (types: real, int, logical, string, vector, matrix).");
      if (kind_of(statement["target"]) != "name")
        this->error(statement, "a constant's name is needed.");
      std::string name = statement["target"]["name"].as_string();
      this->define(statement, name);
      this->declared_constants.push_back(name);
      this->constant_names.insert(name);
      if (statement["operator"].as_string() == "=")
        this->constants.push_back(this->value_of(statement["value"]), name);
    }
  }

  void data_section(const RecipeValue &section)
  {
    const RecipeValue &statements = section["statements"];
    for (size_t i = 0; i < statements.size(); ++i)
    {
      const RecipeValue &statement = statements[i];
      this->refuse_loop(statement, "data");
      if (statement["operator"].as_string().empty())
      {
        // a declaration: the name is a data variable (for calls to other languages)
        if (kind_of(statement["target"]) == "name")
          this->data_names.insert(statement["target"]["name"].as_string());
        continue;
      }
      if (kind_of(statement["target"]) != "name")
        this->error(statement, "a data variable's name is needed.");
      std::string name = statement["target"]["name"].as_string();
      const RecipeValue &value = statement["value"];
      if (this->is_function_value(value))
      {
        // data = cpp::f / r::f(…) (all the data), or y = r::f(…) (one variable)
        bool foreign = this->is_foreign(value);
        if (name != "data" && !foreign)
          this->error(statement, "a C++ function giving the data is written \"data = cpp::f\" (it returns all the data).");
        if (this->data_function)
          this->error(statement, "only one function may give the data.");
        this->data_function = true;
        if (name != "data")
          this->data_names.insert(name);
        RecipeValue entry = RecipeValue::list();
        entry.push_back(foreign ? this->foreign_placeholder(value, "data", name == "data" ? "" : name)
                                : this->function_placeholder(value, "data", "the data"), "data");
        this->data.push_back(entry);
        continue;
      }
      this->data_names.insert(name);
      this->define(statement, name);
      arma::mat numbers;
      if (!this->numeric(value, numbers))
        this->error(value, "data are given by a function (data = cpp::f) or as numbers.");
      this->literal_data.push_back(std::make_pair(name, RecipeValue(numbers)));
    }
  }

  void finish_data()
  {
    if (this->literal_data.empty())
      return;
    if (this->data_function)
      Rcpp::stop(this->file + ": the data are given either by a function (data = cpp::f) or as values, not both.");
    std::string name = "ilike_generated_data_" + std::to_string(++this->generated_count);
    std::string code = "Data " + name + "()\n{\n  Data data;\n";
    for (const auto &item : this->literal_data)
    {
      arma::mat value = item.second.as_matrix();
      code += "  {\n    arma::mat value(" + std::to_string(value.n_rows) + ", " + std::to_string(value.n_cols) + ");\n";
      for (arma::uword i = 0; i < value.n_elem; ++i)
        code += "    value(" + std::to_string(i) + ") = " + number_text(value(i)) + ";\n";
      code += "    data[" + quoted(item.first) + "] = value;\n  }\n";
    }
    code += "  return data;\n}\n\n";
    this->generated_code += code;
    RecipeValue entry = RecipeValue::list();
    entry.push_back(this->placeholder(name, "data"), "data");
    this->data.push_back(entry);
  }

  // the functions of a user-written factor: a bare function (the log
  // density), or density(log = …, simulate = …, gradient = …, second_gradient = …)
  RecipeValue density_functions(const RecipeValue &value, bool prior)
  {
    RecipeValue factor = RecipeValue::list();
    std::string suffix = prior ? "prior" : "likelihood";
    if (this->is_function_value(value))
    {
      factor.push_back(this->function_placeholder(value, "evaluate_log_" + suffix, "the log density"),
                       "evaluate_log_" + suffix);
      return factor;
    }
    if (kind_of(value) != "call" || value["function"].as_string() != "density")
      this->error(value, std::string("a ") + (prior ? "prior" : "factor") +
                  " is a function (cpp::f) or density(log = cpp::f, …).");
    Arguments arguments = this->split_arguments(value);
    std::set<std::string> allowed = {"log", "gradient", "second_gradient"};
    if (prior)
      allowed.insert("simulate");
    this->allow_only(value, arguments, allowed, 0);
    if (arguments.has("log"))
      factor.push_back(this->function_placeholder(arguments.named["log"], "evaluate_log_" + suffix, "log"),
                       "evaluate_log_" + suffix);
    if (arguments.has("simulate"))
      factor.push_back(this->function_placeholder(arguments.named["simulate"], "simulate_prior", "simulate"),
                       "simulate_prior");
    if (arguments.has("gradient"))
      factor.push_back(this->function_placeholder(arguments.named["gradient"], "evaluate_gradient_log_" + suffix, "gradient"), "evaluate_gradient_log_" + suffix);
    if (arguments.has("second_gradient"))
      factor.push_back(this->function_placeholder(arguments.named["second_gradient"], "evaluate_second_gradient_log_" + suffix, "second_gradient"), "evaluate_second_gradient_log_" + suffix);
    if (factor.size() == 0)
      this->error(value, "density(…) needs at least one function.");
    return factor;
  }

  void model_section(const RecipeValue &section)
  {
    const RecipeValue &statements = section["statements"];
    for (size_t i = 0; i < statements.size(); ++i)
    {
      const RecipeValue &statement = statements[i];
      std::string keyword = statement["keyword"].as_string();
      std::string op = statement.contains("operator") ? statement["operator"].as_string() : "";
      if (keyword == "for" || op == "<-" || (op == "~" && keyword.empty() && this->generate))
      {
        // a model statement compiled to generated C++ (I4a)
        this->collect_model_statement(statement, std::vector<ModelLoop>());
        continue;
      }
      const RecipeValue &target = statement["target"];
      if (kind_of(target) == "index")
        this->error(statement, "an indexed target is a model statement (y[i] ~ …) inside a loop.");
      std::string name = target["name"].as_string();
      const RecipeValue &value = statement["value"];

      if (op == "~" && keyword.empty())
      {
        // a built-in prior: theta ~ norm(0, 1)
        if (kind_of(value) != "call" || !role_exists(this->builtin_name(value), BuiltinRole::Prior))
          this->error(value, "after \"~\" a built-in distribution is needed (e.g. theta ~ norm(0, 1)).");
        RecipeValue factor = RecipeValue::list();
        factor.push_back(this->builtin_node(value, true, name), "prior");
        this->add_factor(statement, name, factor);
        continue;
      }
      if (op != "=")
        this->error(statement, "a model statement is \"variable ~ distribution(…)\", \"prior name = …\" or \"factor name = …\".");

      if (keyword.empty() && name == "time")
      {
        // the time index of a state-space model: time = ssm(…)
        if (kind_of(value) != "call" || this->builtin_name(value) != "ssm")
          this->error(value, "time = ssm(index = …, first_index = …, …): the time index of a state-space model.");
        RecipeValue item = RecipeValue::list();
        item.push_back(this->builtin_node(value, false), "ssm");
        this->model_methods.push_back(item);
        continue;
      }
      if (keyword == "prior")
      {
        this->add_factor(statement, name, this->density_functions(value, true));
        continue;
      }
      if (keyword == "transition")
      {
        this->define(statement, "transition " + name);
        this->transitions.push_back(this->transition_model(name, value));
        continue;
      }
      if (keyword == "factor")
      {
        std::string function = (kind_of(value) == "call" && (value["namespace"].as_string().empty() ||
                                                            value["namespace"].as_string() == "ilike"))
                               ? value["function"].as_string() : "";
        if (function == "sbi")
          this->add_factor(statement, name, this->sbi_factor(value));
        else if (function == "linear_gaussian")
          this->add_factor(statement, name, this->gaussian_data_model(value, true));
        else if (function == "nonlinear_gaussian")
          this->add_factor(statement, name, this->gaussian_data_model(value, false));
        else if (!function.empty() && role_exists(function, BuiltinRole::AlgorithmicLikelihood))
          this->add_factor(statement, name, this->algorithmic_likelihood(value));
        else if (!function.empty() && function != "density")
        {
          if (role_exists(function, BuiltinRole::LinearGaussianDataModel))
          {
            RecipeValue factor = RecipeValue::list();
            factor.push_back(this->builtin_node(value, true), "linear_gaussian_data_model");
            this->add_factor(statement, name, factor);
          }
          else
            this->error(value, "unknown kind of factor " + function + " (factors: cpp::f, density(…), sbi(…), "
                        "linear_gaussian(…), nonlinear_gaussian(…), pf(…), kf(…), enkf(…)).");
        }
        else
          this->add_factor(statement, name, this->density_functions(value, false));
        continue;
      }
      this->error(statement, "unknown model statement \"" + keyword + "\" (model statements: variable ~ distribution(…), "
                  "prior name = …, factor name = …, transition variable = …, time = ssm(…)).");
    }
  }

  // --- model statements (I4a) ------------------------------------------------------

  // Whether the model statements need generated code: any loop, deterministic
  // node, indexed target, observed target, or a "~" whose arguments are not
  // all constant (otherwise "theta ~ norm(0, 1)" stays the registry's prior,
  // as before)
  void decide_generation()
  {
    this->generation_decided = true;
    for (size_t t = 0; t < this->trees.size(); ++t)
    {
      const RecipeValue &sections = this->trees[t]["sections"];
      for (size_t s = 0; s < sections.size(); ++s)
      {
        if (sections[s]["section"].as_string() != "model")
          continue;
        const RecipeValue &statements = sections[s]["statements"];
        for (size_t i = 0; i < statements.size(); ++i)
        {
          const RecipeValue &st = statements[i];
          std::string keyword = st["keyword"].as_string();
          std::string op = st.contains("operator") ? st["operator"].as_string() : "";
          if (keyword == "for" || op == "<-")
            this->generate = true;
          if (op == "~" && keyword.empty())
          {
            const RecipeValue &target = st["target"], &value = st["value"];
            if (kind_of(target) == "index" || this->data_names.count(target["name"].as_string()))
              this->generate = true;
            else if (kind_of(value) == "call")
            {
              if (!role_exists(value["function"].as_string(), BuiltinRole::Prior))
                this->generate = true;
              const RecipeValue &arguments = value["arguments"];
              for (size_t k = 0; k < arguments.size(); ++k)
              {
                arma::mat number;
                bool constant = this->numeric(arguments[k], number) ||
                                (kind_of(arguments[k]) == "name" && this->constant_names.count(arguments[k]["name"].as_string()));
                if (!constant)
                  this->generate = true;
              }
            }
          }
        }
      }
    }
  }

  void collect_model_statement(const RecipeValue &statement, std::vector<ModelLoop> loops)
  {
    std::string keyword = statement["keyword"].as_string();
    if (keyword == "for")
    {
      const RecipeValue &range = statement["range"];
      if (kind_of(range) != "binary" || range["operator"].as_string() != ":")
        this->error(statement, "a loop's range is from:to (e.g. for (i in 1:n)).");
      ModelLoop loop;
      loop.index = statement["index"].as_string();
      for (const ModelLoop &l : loops)
        if (l.index == loop.index)
          this->error(statement, "the index " + loop.index + " is already used by an enclosing loop.");
      loop.from = range["left"];
      loop.to = range["right"];
      loops.push_back(loop);
      const RecipeValue &body = statement["body"];
      for (size_t k = 0; k < body.size(); ++k)
        this->collect_model_statement(body[k], loops);
      return;
    }
    std::string op = statement["operator"].as_string();
    if (!keyword.empty() || (op != "~" && op != "<-"))
      this->error(statement, "inside a loop, a statement is \"y[i] ~ distribution(…)\" or \"mu[i] <- expression\".");
    ModelStatement st;
    st.op = op;
    st.node = statement;
    st.loops = loops;
    const RecipeValue &target = statement["target"];
    if (kind_of(target) == "index")
    {
      st.variable = target["target"]["name"].as_string();
      st.indices = target["indices"];
    }
    else
    {
      st.variable = target["name"].as_string();
      st.indices = RecipeValue::list();
    }
    st.value = statement["value"];
    st.observed = this->data_names.count(st.variable) > 0;
    if (this->constant_names.count(st.variable))
      this->error(statement, st.variable + " is a constant, so it cannot be the target of a model statement.");
    if (op == "<-" && st.observed)
      this->error(statement, st.variable + " is data, so it cannot be a deterministic node.");
    if (op == "~")
    {
      if (kind_of(st.value) != "call")
        this->error(st.value, "after \"~\" a distribution is needed (e.g. norm(mu, sigma)).");
      // the factor this statement belongs to, at its first statement's place
      std::string factor_name = st.observed ? st.variable : std::string("model_prior");
      if (!this->model_factor_position.count(factor_name))
      {
        // (named after its data variable, or model_prior: not a new name
        // for the duplicate check, but distinct from the other factors)
        if (this->factor_numbers.count(factor_name))
          this->error(statement, "there is already a factor called " + factor_name + ".");
        this->factor_entries.push_back(RecipeValue::list());
        this->factor_numbers[factor_name] = this->factor_entries.size();
        this->model_factor_position[factor_name] = this->factor_entries.size() - 1;
      }
    }
    this->model_statements.push_back(st);
  }

  // --- generating C++ for the statements

  struct CodeContext
  {
    bool prior = false;                 // in the prior's functions (no runtime data)
    std::set<std::string> indices;      // loop indices in scope
    std::string data = "data";          // how the runtime data are reached ("data" or "(*data)")
  };

  bool is_deterministic(const std::string &name) const
  {
    for (const ModelStatement &st : this->model_statements)
      if (st.op == "<-" && st.variable == name)
        return true;
    return false;
  }

  const ModelStatement &deterministic_statement(const std::string &name) const
  {
    for (const ModelStatement &st : this->model_statements)
      if (st.op == "<-" && st.variable == name)
        return st;
    throw std::runtime_error("internal: no deterministic node " + name);
  }

  // literal data, as a static matrix in the generated code (for priors)
  std::string literal_data_matrix(const std::string &name, const RecipeValue &node)
  {
    for (const auto &item : this->literal_data)
      if (item.first == name)
      {
        std::string symbol = "ilike_data_" + name;
        if (this->generated_code.find("static arma::mat " + symbol + "_make") == std::string::npos)
        {
          arma::mat value = item.second.as_matrix();
          std::string code = "static arma::mat " + symbol + "_make()\n{\n  arma::mat value(" + std::to_string(value.n_rows) +
                             ", " + std::to_string(value.n_cols) + ");\n";
          for (arma::uword k = 0; k < value.n_elem; ++k)
            code += "  value(" + std::to_string(k) + ") = " + number_text(value(k)) + ";\n";
          code += "  return value;\n}\nstatic const arma::mat " + symbol + " = " + symbol + "_make();\n\n";
          this->model_header += code;
        }
        return symbol;
      }
    this->error(node, "the data variable " + name + " is given by a function, so it cannot be used in the prior "
                "(the prior's functions do not receive the data): give it as values in the data section, or as a constant.");
  }

  // a variable as a matrix expression in C++
  std::string variable_matrix(const std::string &name, const RecipeValue &node, const CodeContext &context)
  {
    if (this->constant_names.count(name))
    {
      this->model_constants_used.insert(name);
      return "ilike_constant_" + name;
    }
    if (this->data_names.count(name))
    {
      if (context.prior)
        return this->literal_data_matrix(name, node);
      return context.data + "[" + quoted(name) + "]";
    }
    if (this->is_deterministic(name))
      this->error(node, name + " is a deterministic node: use it with its indices.");
    return "parameters[" + quoted(name) + "]";
  }

  std::string index_code(const RecipeValue &expression, const CodeContext &context)
  {
    return "(int)(" + this->scalar_code(expression, context) + ") - 1";
  }

  // an expression as a double in C++
  std::string scalar_code(const RecipeValue &expression, const CodeContext &context)
  {
    std::string kind = kind_of(expression);
    if (kind == "number")
      return "(" + number_text(expression["value"].as_double()) + ")";
    if (kind == "name")
    {
      if (!expression["namespace"].as_string().empty())
        this->error(expression, "a function cannot be used here.");
      std::string name = expression["name"].as_string();
      if (context.indices.count(name))
        return "(double)" + name;
      if (this->is_deterministic(name))
        return this->deterministic_call(name, RecipeValue::list(), expression, context);
      return "(" + this->variable_matrix(name, expression, context) + ")(0)";
    }
    if (kind == "index")
    {
      const RecipeValue &target = expression["target"];
      if (kind_of(target) != "name")
        this->error(expression, "only a variable can be indexed.");
      std::string name = target["name"].as_string();
      const RecipeValue &indices = expression["indices"];
      if (this->is_deterministic(name))
        return this->deterministic_call(name, indices, expression, context);
      std::string base = this->variable_matrix(name, target, context);
      if (indices.size() == 1)
        return "(" + base + ")(" + this->index_code(indices[0], context) + ")";
      if (indices.size() == 2)
        return "(" + base + ")(" + this->index_code(indices[0], context) + ", " + this->index_code(indices[1], context) + ")";
      this->error(expression, "at most two indices.");
    }
    if (kind == "binary")
    {
      std::string op = expression["operator"].as_string();
      std::string a = this->scalar_code(expression["left"], context), b = this->scalar_code(expression["right"], context);
      if (op == "^")
        return "std::pow(" + a + ", " + b + ")";
      if (op == ":")
        this->error(expression, "a range is only allowed in a loop.");
      return "(" + a + " " + op + " " + b + ")";
    }
    if (kind == "unary")
      return "(-" + this->scalar_code(expression["operand"], context) + ")";
    if (kind == "call")
    {
      static const std::map<std::string, std::string> one = {
        {"exp", "std::exp"}, {"log", "std::log"}, {"sqrt", "std::sqrt"}, {"abs", "std::fabs"},
        {"lgamma", "std::lgamma"}, {"log1p", "std::log1p"}, {"expm1", "std::expm1"},
        {"sin", "std::sin"}, {"cos", "std::cos"}, {"tanh", "std::tanh"}};
      std::string f = expression["function"].as_string();
      const RecipeValue &arguments = expression["arguments"];
      if (!expression["namespace"].as_string().empty())
        this->error(expression, "calls to other languages are not available in model statements (use them in factors).");
      auto argument = [&](size_t k) { return this->scalar_code(arguments[k], context); };
      if (one.count(f))
      {
        if (arguments.size() != 1) this->error(expression, f + " takes one argument.");
        return one.at(f) + "(" + argument(0) + ")";
      }
      if (f == "inv_logit" || f == "ilogit")
      {
        if (arguments.size() != 1) this->error(expression, f + " takes one argument.");
        return "(1.0/(1.0 + std::exp(-" + argument(0) + ")))";
      }
      if (f == "logit")
      {
        if (arguments.size() != 1) this->error(expression, f + " takes one argument.");
        return "std::log(" + argument(0) + "/(1.0 - " + argument(0) + "))";
      }
      if (f == "pow")
      {
        if (arguments.size() != 2) this->error(expression, "pow takes two arguments.");
        return "std::pow(" + argument(0) + ", " + argument(1) + ")";
      }
      if (f == "min" || f == "max")
      {
        if (arguments.size() != 2) this->error(expression, f + " takes two arguments.");
        return "std::" + f + "(" + argument(0) + ", " + argument(1) + ")";
      }
      this->error(expression, "unknown function " + f + " in a model statement (functions: exp, log, sqrt, abs, pow, "
                  "inv_logit, logit, lgamma, log1p, expm1, sin, cos, tanh, min, max).");
    }
    this->error(expression, "this expression cannot be used in a model statement.");
  }

  // a whole vector or matrix (mvnorm's arguments)
  std::string matrix_code(const RecipeValue &expression, const CodeContext &context)
  {
    if (kind_of(expression) == "name" && expression["namespace"].as_string().empty() &&
        !context.indices.count(expression["name"].as_string()))
      return "arma::mat(" + this->variable_matrix(expression["name"].as_string(), expression, context) + ")";
    arma::mat value;
    if (this->numeric(expression, value))
    {
      std::string code = "arma::mat({";
      for (arma::uword r = 0; r < value.n_rows; ++r)
      {
        code += (r ? ", {" : "{");
        for (arma::uword c = 0; c < value.n_cols; ++c)
          code += (c ? ", " : "") + number_text(value(r, c));
        code += "}";
      }
      return code + "})";
    }
    this->error(expression, "a vector or matrix (a variable, a constant or numbers) is needed here.");
  }

  // whether an expression reads data given by a function (directly or through deterministic nodes)
  bool uses_runtime_data(const RecipeValue &expression, std::set<std::string> visited = std::set<std::string>()) const
  {
    std::string kind = kind_of(expression);
    if (kind == "name")
    {
      std::string name = expression["name"].as_string();
      if (this->data_names.count(name))
      {
        for (const auto &item : this->literal_data)
          if (item.first == name)
            return false;
        return true;
      }
      if (this->is_deterministic(name) && !visited.count(name))
      {
        visited.insert(name);
        return this->uses_runtime_data(this->deterministic_statement(name).value, visited);
      }
      return false;
    }
    if (!expression.is_list())
      return false;
    for (size_t k = 0; k < expression.size(); ++k)
      if (expression[k].is_list() && this->uses_runtime_data(expression[k], visited))
        return true;
    return false;
  }

  // the parameters (not data, constants, indices) an expression reads, through deterministic nodes
  void read_parameters(const RecipeValue &expression, std::set<std::string> &found, std::set<std::string> visited = std::set<std::string>()) const
  {
    std::string kind = kind_of(expression);
    if (kind == "name" && expression["namespace"].as_string().empty())
    {
      std::string name = expression["name"].as_string();
      if (this->is_deterministic(name))
      {
        if (!visited.count(name))
        {
          visited.insert(name);
          this->read_parameters(this->deterministic_statement(name).value, found, visited);
        }
      }
      else if (!this->data_names.count(name) && !this->constant_names.count(name))
        found.insert(name);
      return;
    }
    if (!expression.is_list())
      return;
    for (size_t k = 0; k < expression.size(); ++k)
      if (expression[k].is_list())
        this->read_parameters(expression[k], found, visited);
  }

  std::string deterministic_call(const std::string &name, const RecipeValue &indices, const RecipeValue &node, const CodeContext &context)
  {
    const ModelStatement &st = this->deterministic_statement(name);
    if (indices.size() != st.indices.size())
      this->error(node, "the deterministic node " + name + " takes " + std::to_string(st.indices.size()) + " index(es).");
    if (context.prior && this->uses_runtime_data(st.value))
      this->error(node, "the deterministic node " + name + " reads data given by a function, so it cannot be used in the prior.");
    std::string call = "ilike_model_node_" + name + "(parameters, " + (context.prior ? std::string("nullptr") : std::string("&data"));
    for (size_t k = 0; k < indices.size(); ++k)
      call += ", (int)(" + this->scalar_code(indices[k], context) + ")";
    return call + ")";
  }

  // the log density of target ~ dist(arguments), and (for priors) the code
  // that simulates it
  std::string density_code(const ModelStatement &st, const std::string &x, const CodeContext &context)
  {
    const RecipeValue &call = st.value;
    std::string f = call["function"].as_string();
    const RecipeValue &arguments = call["arguments"];
    std::vector<std::string> keys = arguments.names();
    for (const std::string &k : keys)
      if (!k.empty())
        this->error(call, "model statements take the distribution's arguments in order (no names), as in BUGS.");
    static const std::map<std::string, std::pair<size_t, std::string>> densities = {
      {"norm", {2, "dnorm"}}, {"lnorm", {2, "dlnorm"}}, {"gamma", {2, "dgamma"}}, {"unif", {2, "dunif"}},
      {"exp", {1, "dexp"}}, {"pois", {1, "dpois"}}, {"tnorm", {4, "dtnorm"}}, {"beta", {2, "ilike_model_dbeta"}},
      {"bern", {1, "ilike_model_dbern"}}, {"binom", {2, "ilike_model_dbinom"}}};
    if (f == "mvnorm")
    {
      if (arguments.size() != 2)
        this->error(call, "mvnorm takes a mean vector and a covariance matrix.");
      return "dmvnorm(arma::colvec(arma::vectorise(" + x + ")), arma::colvec(arma::vectorise(" + this->matrix_code(arguments[0], context) +
             ")), " + this->matrix_code(arguments[1], context) + ")";
    }
    auto found = densities.find(f);
    if (found == densities.end())
      this->error(call, "unknown distribution " + f + " in a model statement (norm, lnorm, gamma, unif, exp, pois, tnorm, beta, bern, binom, mvnorm).");
    if (arguments.size() != found->second.first)
      this->error(call, f + " takes " + std::to_string(found->second.first) + " argument(s).");
    std::string code = found->second.second + "(" + x;
    for (size_t k = 0; k < arguments.size(); ++k)
      code += ", " + this->scalar_code(arguments[k], context);
    return code + ")";
  }

  std::string random_code(const ModelStatement &st, const CodeContext &context)
  {
    const RecipeValue &call = st.value;
    std::string f = call["function"].as_string();
    const RecipeValue &arguments = call["arguments"];
    static const std::map<std::string, std::string> samplers = {
      {"norm", "rnorm"}, {"lnorm", "rlnorm"}, {"gamma", "rgamma"}, {"unif", "runif"}, {"exp", "rexp"},
      {"pois", "(double)rpois"}, {"tnorm", "rtnorm"}, {"beta", "ilike_model_rbeta"}};
    if (f == "mvnorm")
      return "arma::mat(rmvnorm(rng, arma::colvec(arma::vectorise(" + this->matrix_code(arguments[0], context) + ")), " +
             this->matrix_code(arguments[1], context) + "))";
    auto found = samplers.find(f);
    if (found == samplers.end())
      this->error(call, "the parameter " + st.variable + " has the distribution " + f +
                  ", which cannot be simulated in a prior (a discrete parameter needs element moves, I4c).");
    std::string code = found->second + "(rng";
    for (size_t k = 0; k < arguments.size(); ++k)
      code += ", " + this->scalar_code(arguments[k], context);
    return code + ")";
  }

  // the target of a statement as C++ (an element, a scalar, or a whole vector for mvnorm)
  std::string target_code(const ModelStatement &st, const CodeContext &context, bool writable)
  {
    std::string base = st.observed ? (context.data + "[" + quoted(st.variable) + "]") : ("parameters[" + quoted(st.variable) + "]");
    if (st.value["function"].as_string() == "mvnorm")
    {
      if (st.indices.size() != 0)
        this->error(st.node, "mvnorm's target is a whole vector (no index).");
      return base;
    }
    if (st.indices.size() == 0)
      return "(" + base + ")(0)";
    if (st.indices.size() == 1)
      return "(" + base + ")(" + this->index_code(st.indices[0], context) + ")";
    if (st.indices.size() == 2)
      return "(" + base + ")(" + this->index_code(st.indices[0], context) + ", " + this->index_code(st.indices[1], context) + ")";
    this->error(st.node, "at most two indices.");
  }

  // loops around a body: for (int i = from; i <= to; ++i) { … }
  std::string loops_code(const ModelStatement &st, CodeContext &context, const std::string &body, const std::string &indent)
  {
    std::string open, close, pad = indent;
    for (const ModelLoop &loop : st.loops)
    {
      std::string from = this->scalar_code(loop.from, context), to = this->scalar_code(loop.to, context);
      open += pad + "for (int " + loop.index + " = (int)" + from + "; " + loop.index + " <= (int)" + to + "; ++" + loop.index + ")\n" + pad + "{\n";
      context.indices.insert(loop.index);
      close = pad + "}\n" + close;
      pad += "  ";
    }
    std::string inner;
    std::stringstream lines(body);
    std::string line;
    while (std::getline(lines, line))
      inner += pad + line + "\n";
    return open + inner + close;
  }

  void finish_model_statements()
  {
    if (this->model_statements.empty())
      return;
    std::string code = R"(// model statements (generated by the ilike language, I4a)
inline double ilike_model_dbeta(double x, double a, double b)
{
  if (x <= 0.0 || x >= 1.0) return -std::numeric_limits<double>::infinity();
  return std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b) + (a - 1.0)*std::log(x) + (b - 1.0)*std::log(1.0 - x);
}
inline double ilike_model_rbeta(RandomNumberGenerator &rng, double a, double b)
{
  double x = rgamma(rng, a, 1.0), y = rgamma(rng, b, 1.0);
  return x/(x + y);
}
inline double ilike_model_dbern(double x, double p)
{
  if (x == 1.0) return std::log(p);
  if (x == 0.0) return std::log(1.0 - p);
  return -std::numeric_limits<double>::infinity();
}
inline double ilike_model_dbinom(double x, double n, double p)
{
  if (x < 0.0 || x > n || x != std::floor(x)) return -std::numeric_limits<double>::infinity();
  return std::lgamma(n + 1.0) - std::lgamma(x + 1.0) - std::lgamma(n - x + 1.0) + x*std::log(p) + (n - x)*std::log(1.0 - p);
}

)";
    // deterministic nodes: one function each (data may be absent, in the prior)
    std::set<std::string> deterministic_done;
    std::string node_code;
    for (const ModelStatement &st : this->model_statements)
    {
      if (st.op != "<-")
        continue;
      if (deterministic_done.count(st.variable))
        this->error(st.node, "the deterministic node " + st.variable + " is defined more than once.");
      deterministic_done.insert(st.variable);
      CodeContext context;
      context.data = "(*data)";
      std::string signature = "inline double ilike_model_node_" + st.variable + "(const Parameters &parameters, const Data *data";
      for (size_t k = 0; k < st.indices.size(); ++k)
      {
        if (kind_of(st.indices[k]) != "name")
          this->error(st.node, "a deterministic node's indices are the loop's indices (e.g. mu[i] <- …).");
        std::string index = st.indices[k]["name"].as_string();
        signature += ", int " + index;
        context.indices.insert(index);
      }
      node_code += signature + ");\n";
      std::string body = signature + ")\n{\n  return " + this->scalar_code(st.value, context) + ";\n}\n\n";
      this->deferred_node_code += body;
    }
    code += node_code + "\n";

    // the likelihood factors: one per observed variable
    std::vector<std::string> observed_done;
    for (const ModelStatement &first : this->model_statements)
    {
      if (first.op != "~" || !first.observed ||
          std::find(observed_done.begin(), observed_done.end(), first.variable) != observed_done.end())
        continue;
      observed_done.push_back(first.variable);
      std::string function = "ilike_model_log_likelihood_" + first.variable;
      std::string body = "double log_density = 0.0;\n";
      for (const ModelStatement &st : this->model_statements)
      {
        if (st.op != "~" || st.variable != first.variable)
          continue;
        // (the loop indices in scope for the target and the arguments)
        CodeContext inner_context;
        for (const ModelLoop &loop : st.loops) inner_context.indices.insert(loop.index);
        std::string statement_body = "log_density += " + this->density_code(st, this->target_code(st, inner_context, false), inner_context) + ";\n";
        CodeContext fresh;
        body += this->loops_code(st, fresh, statement_body, "");
      }
      body += "return log_density;\n";
      std::string indented;
      std::stringstream lines(body);
      std::string line;
      while (std::getline(lines, line)) indented += "  " + line + "\n";
      this->deferred_function_code += "double " + function + "(const Parameters &parameters, const Data &data)\n{\n" + indented + "}\n\n";
      RecipeValue factor = RecipeValue::list();
      factor.push_back(this->placeholder(function, "evaluate_log_likelihood"), "evaluate_log_likelihood");
      this->factor_entries[this->model_factor_position.at(first.variable)] = factor;
    }

    // the joint prior: its log density, and its simulation in dependency order
    std::vector<size_t> prior;
    for (size_t k = 0; k < this->model_statements.size(); ++k)
      if (this->model_statements[k].op == "~" && !this->model_statements[k].observed)
        prior.push_back(k);
    if (!prior.empty())
    {
      std::set<std::string> prior_variables;
      for (size_t k : prior)
        prior_variables.insert(this->model_statements[k].variable);
      // the parameters each statement reads (arguments, indices and loop bounds)
      std::map<size_t, std::set<std::string>> reads;
      for (size_t k : prior)
      {
        const ModelStatement &st = this->model_statements[k];
        std::set<std::string> found;
        this->read_parameters(st.value, found);
        for (size_t i = 0; i < st.indices.size(); ++i) this->read_parameters(st.indices[i], found);
        for (const ModelLoop &loop : st.loops) { this->read_parameters(loop.from, found); this->read_parameters(loop.to, found); }
        for (const ModelLoop &loop : st.loops) found.erase(loop.index);
        for (const std::string &name : found)
          if (!prior_variables.count(name))
            this->error(st.node, "the prior of " + st.variable + " depends on " + name +
                        ", which has no distribution in the model statements (and is not data or a constant).");
        found.erase(st.variable); // (an element may depend on earlier elements, e.g. x[t] on x[t-1])
        reads[k] = found;
      }
      // order: a statement after every statement of the variables it reads;
      // statements of one variable in the file's order
      std::vector<size_t> order;
      std::set<size_t> placed;
      while (order.size() < prior.size())
      {
        bool progress = false;
        for (size_t k : prior)
        {
          if (placed.count(k))
            continue;
          bool ready = true;
          for (size_t j : prior)
            if (!placed.count(j) && j != k &&
                (reads[k].count(this->model_statements[j].variable) ||
                 (this->model_statements[j].variable == this->model_statements[k].variable && j < k)))
              ready = false;
          if (ready)
          {
            order.push_back(k);
            placed.insert(k);
            progress = true;
          }
        }
        if (!progress)
          this->error(this->model_statements[prior[0]].node, "the prior's statements depend on each other in a cycle.");
      }

      std::string log_body = "double log_density = 0.0;\n";
      std::string simulate_body = "Parameters parameters;\n";
      std::set<std::string> allocated;
      for (size_t k : order)
      {
        const ModelStatement &st = this->model_statements[k];
        CodeContext loop_context;
        loop_context.prior = true;
        for (const ModelLoop &loop : st.loops) loop_context.indices.insert(loop.index);
        std::string x = this->target_code(st, loop_context, false);
        CodeContext fresh;
        fresh.prior = true;
        log_body += this->loops_code(st, fresh, "log_density += " + this->density_code(st, x, loop_context) + ";\n", "");

        if (!allocated.count(st.variable))
        {
          allocated.insert(st.variable);
          // the variable's size: the largest index of its statements
          if (st.value["function"].as_string() != "mvnorm")
          {
            std::string rows = "1", cols = "1";
            for (const ModelStatement &other : this->model_statements)
            {
              if (other.op != "~" || other.variable != st.variable)
                continue;
              CodeContext size_context;
              size_context.prior = true;
              for (size_t i = 0; i < other.indices.size() && i < 2; ++i)
              {
                std::string extent;
                const RecipeValue &index = other.indices[i];
                if (kind_of(index) == "name")
                {
                  std::string name = index["name"].as_string();
                  for (const ModelLoop &loop : other.loops)
                    if (loop.index == name)
                      extent = "(int)" + this->scalar_code(loop.to, size_context);
                }
                if (extent.empty())
                {
                  arma::mat number;
                  if (!this->numeric(index, number))
                    this->error(other.node, "a parameter's index in its prior is a loop index or a number (e.g. x[1], x[t]).");
                  extent = std::to_string((int)number(0));
                }
                std::string &dimension = (i == 0) ? rows : cols;
                dimension = "std::max(" + dimension + ", " + extent + ")";
              }
            }
            simulate_body += "parameters[" + quoted(st.variable) + "] = arma::mat(" + rows + ", " + cols + ", arma::fill::zeros);\n";
          }
        }
        std::string assign;
        if (st.value["function"].as_string() == "mvnorm")
          assign = "parameters[" + quoted(st.variable) + "] = " + this->random_code(st, loop_context) + ";\n";
        else
          assign = x + " = " + this->random_code(st, loop_context) + ";\n";
        CodeContext fresh2;
        fresh2.prior = true;
        simulate_body += this->loops_code(st, fresh2, assign, "");
      }
      log_body += "return log_density;\n";
      simulate_body += "return parameters;\n";
      auto indent = [](const std::string &body) {
        std::string out, line;
        std::stringstream lines(body);
        while (std::getline(lines, line)) out += "  " + line + "\n";
        return out;
      };
      this->deferred_function_code += "double ilike_model_log_prior(const Parameters &parameters)\n{\n" + indent(log_body) + "}\n\n";
      this->deferred_function_code += "Parameters ilike_model_simulate_prior(RandomNumberGenerator &rng)\n{\n" + indent(simulate_body) + "}\n\n";
      RecipeValue factor = RecipeValue::list();
      factor.push_back(this->placeholder("ilike_model_log_prior", "evaluate_log_prior"), "evaluate_log_prior");
      factor.push_back(this->placeholder("ilike_model_simulate_prior", "simulate_prior"), "simulate_prior");
      this->factor_entries[this->model_factor_position.at("model_prior")] = factor;
    }

    // a deterministic node defined but never used is harmless; a parameter
    // read by a likelihood but given no distribution may be given by another factor
    this->model_code = code;
  }

  // a matrix argument of a Gaussian model: a function (cpp::f, r::f, …), or
  // numbers, for which a function is generated (taking the parameters if
  // with_parameters)
  RecipeValue matrix_placeholder(const RecipeValue &expression, const std::string &what, bool with_parameters,
                                 const std::string &role)
  {
    if (this->is_function_value(expression))
      return this->function_placeholder(expression, role, what);
    arma::mat value;
    if (!this->numeric(expression, value))
      this->error(expression, what + " is a matrix (numbers) or a function (cpp::f).");
    std::string name = "ilike_generated_matrix_" + std::to_string(++this->generated_count);
    std::string code = "arma::mat " + name + "(" + (with_parameters ? "const Parameters &parameters" : "") + ")\n{\n" +
                       "  arma::mat value(" + std::to_string(value.n_rows) + ", " + std::to_string(value.n_cols) + ");\n";
    for (arma::uword k = 0; k < value.n_elem; ++k)
      code += "  value(" + std::to_string(k) + ") = " + number_text(value(k)) + ";\n";
    code += "  return value;\n}\n\n";
    this->generated_code += code;
    return this->placeholder(name, role);
  }

  // TRUE if a function (cpp::f) among these takes no arguments: the
  // generated matrix functions then take none either, so that the recipe's
  // "type" (with or without parameters) is the same for all of them
  bool functions_take_parameters(const std::vector<RecipeValue> &expressions)
  {
    for (const RecipeValue &e : expressions)
      if (this->is_function_reference(e))
      {
        std::string signature = find_cpp_function_signature(this->user_code, e["name"].as_string());
        size_t open = signature.find('('), close = signature.rfind(')');
        if (open != std::string::npos && close != std::string::npos)
        {
          std::string inside = signature.substr(open + 1, close - open - 1);
          if (inside.find_first_not_of(" \t\r\n") == std::string::npos)
            return false;
        }
      }
    return true;
  }

  // linear_gaussian(data = y, state = x, matrix = …, covariance = …) or
  // nonlinear_gaussian(data = y, function = cpp::h, covariance = …)
  RecipeValue gaussian_data_model(const RecipeValue &call, bool linear)
  {
    Arguments arguments = this->split_arguments(call);
    std::vector<std::string> order = linear ? std::vector<std::string>{"data", "state", "matrix", "covariance"}
                                            : std::vector<std::string>{"data", "function", "covariance"};
    this->allow_only(call, arguments, std::set<std::string>(order.begin(), order.end()), order.size());
    for (size_t k = 0; k < arguments.positional.size(); ++k)
    {
      if (arguments.has(order[k]))
        this->error(call, call["function"].as_string() + ": " + order[k] + " is given twice.");
      arguments.named[order[k]] = arguments.positional[k];
    }
    for (const std::string &a : order)
      if (!arguments.has(a))
        this->error(call, call["function"].as_string() + ": " + a + " is missing.");
    auto variable_name = [&](const std::string &a) {
      const RecipeValue &e = arguments.named[a];
      if (kind_of(e) == "name" && !this->is_function_reference(e))
        return e["name"].as_string();
      if (kind_of(e) == "vector")
      {
        std::string names;
        for (size_t k = 0; k < e["elements"].size(); ++k)
        {
          if (kind_of(e["elements"][k]) != "name")
            this->error(e, call["function"].as_string() + ": " + a + " is a variable name (or a list of names).");
          names = names + (k ? ";" : "") + e["elements"][k]["name"].as_string();
        }
        return names;
      }
      this->error(e, call["function"].as_string() + ": " + a + " is a variable name.");
    };
    std::string data = variable_name("data");
    RecipeValue factor = RecipeValue::list();
    if (linear)
    {
      std::string state = variable_name("state");
      arma::mat a, b;
      if (this->numeric(arguments.named["matrix"], a) && this->numeric(arguments.named["covariance"], b))
      {
        // the one-line built-in
        RecipeValue node = RecipeValue::list();
        node.push_back(RecipeValue("linear_gaussian_data_model"), "type");
        node.push_back(RecipeValue(data), "variables");
        RecipeValue parameters = RecipeValue::list();
        parameters.push_back(RecipeValue(state));
        parameters.push_back(this->value_of(arguments.named["matrix"]));
        parameters.push_back(this->value_of(arguments.named["covariance"]));
        node.push_back(parameters, "parameters");
        factor.push_back(node, "linear_gaussian_data_model");
        return factor;
      }
      bool with_parameters = this->functions_take_parameters({arguments.named["matrix"], arguments.named["covariance"]});
      factor.push_back(this->placeholder(this->generate_string_function(data), "linear_gaussian_data_variable"),
                       "linear_gaussian_data_variable");
      factor.push_back(this->placeholder(this->generate_string_function(state), "linear_gaussian_data_state_variable"),
                       "linear_gaussian_data_state_variable");
      factor.push_back(this->matrix_placeholder(arguments.named["matrix"], "matrix", with_parameters, "linear_gaussian_data_matrix"), "linear_gaussian_data_matrix");
      factor.push_back(this->matrix_placeholder(arguments.named["covariance"], "covariance", with_parameters, "linear_gaussian_data_covariance"), "linear_gaussian_data_covariance");
      return factor;
    }
    bool with_parameters = this->functions_take_parameters({arguments.named["covariance"]});
    factor.push_back(this->placeholder(this->generate_string_function(data), "nonlinear_gaussian_data_variable"),
                     "nonlinear_gaussian_data_variable");
    factor.push_back(this->function_placeholder(arguments.named["function"], "nonlinear_gaussian_data_function", "function"), "nonlinear_gaussian_data_function");
    factor.push_back(this->matrix_placeholder(arguments.named["covariance"], "covariance", with_parameters, "nonlinear_gaussian_data_covariance"), "nonlinear_gaussian_data_covariance");
    return factor;
  }

  // transition x = linear_gaussian(matrix = …, covariance = …),
  // nonlinear_gaussian(function = cpp::f, covariance = …),
  // density(simulate = cpp::f, log = cpp::g), or a built-in (norm_rw(x, …))
  RecipeValue transition_model(const std::string &variable, const RecipeValue &value)
  {
    RecipeValue entry = RecipeValue::list();
    if (kind_of(value) != "call")
      this->error(value, "a transition is linear_gaussian(…), nonlinear_gaussian(…), density(…) or a built-in (e.g. norm_rw(x, 1)).");
    std::string function = value["function"].as_string();
    Arguments arguments = this->split_arguments(value);
    if (function == "linear_gaussian" || function == "nonlinear_gaussian")
    {
      bool linear = (function == "linear_gaussian");
      std::vector<std::string> order = linear ? std::vector<std::string>{"matrix", "covariance"}
                                              : std::vector<std::string>{"function", "covariance"};
      this->allow_only(value, arguments, std::set<std::string>(order.begin(), order.end()), order.size());
      for (size_t k = 0; k < arguments.positional.size(); ++k)
        arguments.named[order[k]] = arguments.positional[k];
      for (const std::string &a : order)
        if (!arguments.has(a))
          this->error(value, function + ": " + a + " is missing.");
      std::string prefix = linear ? "linear_gaussian_transition_" : "nonlinear_gaussian_transition_";
      arma::mat a, b;
      if (linear && this->numeric(arguments.named["matrix"], a) && this->numeric(arguments.named["covariance"], b))
      {
        RecipeValue node = RecipeValue::list();
        node.push_back(RecipeValue("linear_gaussian_transition_model"), "type");
        node.push_back(RecipeValue(variable), "variables");
        RecipeValue parameters = RecipeValue::list();
        parameters.push_back(this->value_of(arguments.named["matrix"]));
        parameters.push_back(this->value_of(arguments.named["covariance"]));
        node.push_back(parameters, "parameters");
        entry.push_back(node, "linear_gaussian_transition_model");
        return entry;
      }
      bool with_parameters = linear ? this->functions_take_parameters({arguments.named["matrix"], arguments.named["covariance"]})
                                    : this->functions_take_parameters({arguments.named["covariance"]});
      entry.push_back(this->placeholder(this->generate_string_function(variable), prefix + "variable"), prefix + "variable");
      if (linear)
        entry.push_back(this->matrix_placeholder(arguments.named["matrix"], "matrix", with_parameters, prefix + "matrix"), prefix + "matrix");
      else
        entry.push_back(this->function_placeholder(arguments.named["function"], prefix + "function", "function"),
                        prefix + "function");
      entry.push_back(this->matrix_placeholder(arguments.named["covariance"], "covariance", with_parameters, prefix + "covariance"), prefix + "covariance");
      return entry;
    }
    if (function == "density")
    {
      this->allow_only(value, arguments, {"simulate", "log"}, 0);
      if (!arguments.has("simulate"))
        this->error(value, "a transition written by the user needs density(simulate = cpp::f, log = cpp::g) (log is optional for some algorithms).");
      entry.push_back(this->function_placeholder(arguments.named["simulate"], "simulate_transition_model", "simulate"),
                      "simulate_transition_model");
      if (arguments.has("log"))
        entry.push_back(this->function_placeholder(arguments.named["log"], "evaluate_log_transition_model", "log"),
                        "evaluate_log_transition_model");
      return entry;
    }
    std::string type = this->builtin_name(value);
    if (!role_exists(type, BuiltinRole::TransitionModel))
      find_builtin(type, BuiltinRole::TransitionModel);
    // the variable may be repeated as the first argument (norm_rw(x, 1)) or left out (norm_rw(1))
    bool repeated = !arguments.positional.empty() && kind_of(arguments.positional[0]) == "name" &&
                    !this->is_function_reference(arguments.positional[0]) &&
                    arguments.positional[0]["name"].as_string() == variable;
    entry.push_back(repeated ? this->builtin_node(value, true) : this->builtin_node(value, true, variable), "transition_model");
    return entry;
  }

  // pf(inner, particles = …), kf(inner), enkf(inner, …): a filter on an
  // imported file (the binding compiles it; settings missing here are taken
  // from the inner file's algorithm section, dev/ilike-language.md 4.7)
  RecipeValue algorithmic_likelihood(const RecipeValue &call)
  {
    std::string type = this->builtin_name(call);
    const RecipeValue &arguments = call["arguments"];
    std::vector<std::string> keys = arguments.names();
    RecipeValue module;
    std::string algorithm;
    RecipeValue parameters = RecipeValue::list();
    for (size_t k = 0; k < arguments.size(); ++k)
    {
      std::string key = (k < keys.size()) ? keys[k] : "";
      if ((key.empty() && module.is_null()) || key == "model")
      {
        module = arguments[k];
        continue;
      }
      if (key == "algorithm")
      {
        if (kind_of(arguments[k]) != "string" && kind_of(arguments[k]) != "name")
          this->error(arguments[k], type + ": algorithm is the name of the inner file's algorithm.");
        algorithm = kind_of(arguments[k]) == "string" ? arguments[k]["value"].as_string() : arguments[k]["name"].as_string();
        continue;
      }
      parameters.push_back(this->value_of(arguments[k]), key);
    }
    if (module.is_null() || kind_of(module) != "name" || this->is_function_reference(module))
      this->error(call, type + ": the first argument is an imported file (import \"inner.ilike\" as inner).");
    std::string alias = module["name"].as_string();
    auto found = this->modules.find(alias);
    if (found == this->modules.end())
      this->error(module, type + ": " + alias + " is not imported (import \"file.ilike\" as " + alias + ").");
    RecipeValue placeholder = RecipeValue::list();
    placeholder.push_back(RecipeValue(alias), "ilike_module");
    placeholder.push_back(RecipeValue(found->second.first), "path");
    placeholder.push_back(RecipeValue(found->second.second), "importer");
    placeholder.push_back(RecipeValue(algorithm), "algorithm");
    RecipeValue node = RecipeValue::list();
    node.push_back(RecipeValue(type), "type");
    node.push_back(RecipeValue(found->second.first), "model_file");
    node.push_back(placeholder, "model");
    node.push_back(parameters, "parameters");
    RecipeValue factor = RecipeValue::list();
    factor.push_back(node, "algorithmic_likelihood");
    return factor;
  }

  // sbi(method(…), simulate = cpp::f, data = y, summaries = cpp::g)
  RecipeValue sbi_factor(const RecipeValue &call)
  {
    Arguments arguments = this->split_arguments(call);
    this->allow_only(call, arguments, {"method", "simulate", "data", "summaries"}, 1);
    RecipeValue method = arguments.positional.empty() ? (arguments.has("method") ? arguments.named["method"] : RecipeValue())
                                                      : arguments.positional[0];
    if (method.is_null())
      this->error(call, "sbi: the method is missing (e.g. sbi(sl(100), …)).");
    if (!role_exists(this->builtin_name(method), BuiltinRole::SBILikelihood))
      find_builtin(this->builtin_name(method), BuiltinRole::SBILikelihood); // (the registry's message)
    if (!arguments.has("simulate"))
      this->error(call, "sbi: the simulator is missing (simulate = cpp::f).");
    if (!arguments.has("data"))
      this->error(call, "sbi: the data variable is missing (data = y).");
    RecipeValue factor = RecipeValue::list();
    factor.push_back(this->function_placeholder(arguments.named["simulate"], "simulate_data_model", "simulate"),
                     "simulate_data_model");
    factor.push_back(this->builtin_node(method, false), "sbi_likelihood");
    const RecipeValue &data_variable = arguments.named["data"];
    std::string variables;
    if (kind_of(data_variable) == "name")
      variables = data_variable["name"].as_string();
    else if (kind_of(data_variable) == "vector")
    {
      const RecipeValue &elements = data_variable["elements"];
      for (size_t i = 0; i < elements.size(); ++i)
        variables = variables + (i ? ";" : "") + elements[i]["name"].as_string();
    }
    else
      this->error(data_variable, "sbi: data is the name of the data variable (or a list of names).");
    factor.push_back(this->placeholder(this->generate_string_function(variables), "data_variable"), "data_variable");
    if (arguments.has("summaries"))
      factor.push_back(this->function_placeholder(arguments.named["summaries"], "summary_statistics", "summaries"),
                       "summary_statistics");
    return factor;
  }

  // --- algorithms ----------------------------------------------------------------------

  struct Move
  {
    std::string name;
    std::string prefix;  // m, mh, imh, unadjusted
    int code;            // in order_of_mcmc: mh 1, imh 2, m 3, unadjusted 4
    RecipeValue entry;
  };

  Move lower_move(const RecipeValue &statement, const std::string &name, const RecipeValue &call)
  {
    static const std::map<std::string, std::pair<std::string, int>> kinds = {
      {"metropolis", {"m", 3}}, {"metropolis_hastings", {"mh", 1}},
      {"independent_mh", {"imh", 2}}, {"unadjusted", {"unadjusted", 4}}};
    if (kind_of(call) != "call" || !kinds.count(call["function"].as_string()))
      this->error(statement["value"], "a move is metropolis(…), metropolis_hastings(…), independent_mh(…) or unadjusted(…).");
    std::string kind = call["function"].as_string();
    Move move;
    move.name = name;
    move.prefix = kinds.at(kind).first;
    move.code = kinds.at(kind).second;
    const std::string &p = move.prefix;

    Arguments arguments = this->split_arguments(call);
    this->allow_only(call, arguments, {"proposal", "log", "factors", "cut_factors", "adapt", "smc_adapt"}, 1);
    RecipeValue proposal = !arguments.positional.empty() ? arguments.positional[0]
                         : (arguments.has("proposal") ? arguments.named["proposal"] : RecipeValue());
    if (proposal.is_null())
      this->error(call, kind + ": the proposal is missing.");

    RecipeValue entry = RecipeValue::list();
    std::string variable;
    BuiltinRole role = (p == "m") ? BuiltinRole::MProposal : (p == "mh") ? BuiltinRole::MHProposal
                     : (p == "imh") ? BuiltinRole::IMHProposal : BuiltinRole::UnadjustedProposal;
    if (kind_of(proposal) == "call" && proposal["function"].as_string() != "density" && !this->is_function_value(proposal))
    {
      std::string type = this->builtin_name(proposal);
      if (!role_exists(type, role))
        // the registry's message: what it is, and where it may be used
        find_builtin(type, role);
      RecipeValue node = this->builtin_node(proposal, true);
      variable = node["variables"].as_string();
      entry.push_back(node, p + "_proposal");
      if (arguments.has("log"))
        this->error(arguments.named["log"], kind + ": log goes with a proposal written by the user.");
    }
    else
    {
      RecipeValue simulate, log;
      if (this->is_function_value(proposal))
      {
        simulate = proposal;
        if (arguments.has("log")) log = arguments.named["log"];
      }
      else if (kind_of(proposal) == "call")
      {
        Arguments density = this->split_arguments(proposal);
        this->allow_only(proposal, density, {"simulate", "log"}, 0);
        if (density.has("simulate")) simulate = density.named["simulate"];
        if (density.has("log")) log = density.named["log"];
      }
      if (simulate.is_null())
        this->error(proposal, kind + ": the proposal is a built-in (e.g. norm_rw(theta, 1)) or a function (cpp::f).");
      entry.push_back(this->function_placeholder(simulate, "simulate_" + p + "_proposal", "the proposal"),
                      "simulate_" + p + "_proposal");
      bool needs_log = (p == "mh" || p == "imh");
      if (needs_log && log.is_null())
        this->error(call, kind + ": a proposal written by the user needs its log density (log = cpp::g).");
      if (!needs_log && !log.is_null())
        this->error(log, kind + ": no log density is needed (the proposal is " + (p == "m" ? "symmetric" : "not corrected") + ").");
      if (!log.is_null())
        entry.push_back(this->function_placeholder(log, "evaluate_log_" + p + "_proposal", "log"),
                        "evaluate_log_" + p + "_proposal");
    }

    if (arguments.has("factors"))
      entry.push_back(this->placeholder(this->generate_factor_index(arguments.named["factors"]), p + "_factor_index"),
                      p + "_factor_index");
    if (arguments.has("cut_factors"))
      entry.push_back(this->placeholder(this->generate_factor_index(arguments.named["cut_factors"]), p + "_cut_factor_index"),
                      p + "_cut_factor_index");
    if (arguments.has("adapt"))
    {
      if (p != "m" && p != "mh")
        this->error(arguments.named["adapt"], kind + ": adaptation is available for metropolis and metropolis_hastings moves.");
      this->lower_adaptor(arguments.named["adapt"], p, variable, entry);
    }
    if (arguments.has("smc_adapt"))
    {
      const RecipeValue &adaptor = arguments.named["smc_adapt"];
      if (variable.empty())
        this->error(adaptor, kind + ": smc_adapt needs a built-in proposal (it sets its covariance).");
      entry.push_back(this->builtin_node(adaptor, true, variable), "smc_adaptor");
    }
    move.entry = entry;
    return move;
  }

  // adapt = scaled_am(…) (its variable defaults to the proposal's),
  // acceptance_rate(parameter, …), cpp::f, or averaged(cpp::f)
  void lower_adaptor(const RecipeValue &adaptor, const std::string &p, const std::string &variable, RecipeValue &entry)
  {
    if (this->is_function_reference(adaptor))
    {
      entry.push_back(this->function_placeholder(adaptor, "adapt_" + p + "_proposal", "the adaptation"),
                      "adapt_" + p + "_proposal");
      return;
    }
    if (kind_of(adaptor) != "call")
      this->error(adaptor, "adapt is a built-in adaptor (e.g. scaled_am()) or a function (cpp::f).");
    std::string type = this->builtin_name(adaptor);
    if (type == "averaged" && adaptor["arguments"].size() == 0)
    {
      // averaged() alone: refused by the engine, which says it goes with a
      // user's adaptation function (as the old format)
      RecipeValue node = RecipeValue::list();
      node.push_back(RecipeValue("averaged"), "type");
      node.push_back(RecipeValue(""), "variables");
      node.push_back(RecipeValue::list(), "parameters");
      entry.push_back(node, p + "_adaptor");
      return;
    }
    if (type == "averaged")
    {
      Arguments arguments = this->split_arguments(adaptor);
      if (arguments.positional.size() != 1 || !arguments.named.empty() || !this->is_function_reference(arguments.positional[0]))
        this->error(adaptor, "averaged(cpp::f): averaging of a user's adaptation function.");
      entry.push_back(this->function_placeholder(arguments.positional[0], "adapt_" + p + "_proposal", "the adaptation"),
                      "adapt_" + p + "_proposal");
      RecipeValue node = RecipeValue::list();
      node.push_back(RecipeValue("averaged"), "type");
      node.push_back(RecipeValue(""), "variables");
      node.push_back(RecipeValue::list(), "parameters");
      entry.push_back(node, p + "_adaptor");
      return;
    }
    const Builtin &builtin = find_builtin(type, BuiltinRole::MCMCAdaptor);
    bool first_is_variable = !builtin.arguments.empty() && builtin.arguments[0].name == "variable";
    Arguments arguments = this->split_arguments(adaptor);
    if (first_is_variable && !arguments.has("variable"))
    {
      if (variable.empty())
        this->error(adaptor, type + ": give the variable (variable = …), since the proposal is written by the user.");
      entry.push_back(this->builtin_node(adaptor, true, variable), p + "_adaptor");
    }
    else
      entry.push_back(this->builtin_node(adaptor, true), p + "_adaptor");
  }

  void algorithm_section(const RecipeValue &section)
  {
    std::string name = section["name"].as_string();
    if (name.empty())
      name = "main";
    if (this->algorithms.contains(name))
      this->error(section, "algorithm " + name + " is defined more than once.");

    std::vector<Move> moves;
    std::vector<std::string> sweep;
    bool has_sweep = false;
    RecipeValue methods = RecipeValue::list();
    RecipeValue is_proposals = RecipeValue::list();
    RecipeValue transition_proposals = RecipeValue::list();
    RecipeValue settings = RecipeValue::list();

    auto method = [&](const std::string &key, const RecipeValue &node) {
      RecipeValue item = RecipeValue::list();
      item.push_back(node, key);
      methods.push_back(item);
    };

    const RecipeValue &statements = section["statements"];
    for (size_t i = 0; i < statements.size(); ++i)
    {
      const RecipeValue &statement = statements[i];
      this->refuse_loop(statement, "an algorithm section");
      std::string keyword = statement["keyword"].as_string();
      std::string target = statement["target"]["name"].as_string();
      const RecipeValue &value = statement["value"];
      if (statement["operator"].as_string() != "=")
        this->error(statement, "an algorithm holds settings (name = …) and moves (move name = …).");

      if (keyword == "move" || (keyword.empty() && target == "move"))
      {
        std::string move_name = (keyword == "move") ? target : "move" + std::to_string(moves.size() + 1);
        for (const Move &m : moves)
          if (m.name == move_name)
            this->error(statement, "move " + move_name + " is defined more than once.");
        moves.push_back(this->lower_move(statement, move_name, value));
        continue;
      }
      if (!keyword.empty())
        this->error(statement, "unknown kind of statement \"" + keyword + "\" in an algorithm (moves are \"move name = …\").");

      if (target == "sweep")
      {
        if (kind_of(value) != "call" || value["function"].as_string() != "deterministic")
          this->error(value, "sweep = deterministic(m1, m2, …): the moves in order.");
        const RecipeValue &arguments = value["arguments"];
        for (size_t k = 0; k < arguments.size(); ++k)
        {
          if (kind_of(arguments[k]) != "name")
            this->error(arguments[k], "sweep: the names of moves are needed.");
          sweep.push_back(arguments[k]["name"].as_string());
        }
        has_sweep = true;
      }
      else if (target == "termination")
      {
        std::string type = this->builtin_name(value);
        if (role_exists(type, BuiltinRole::MCMCTermination))
          method("mcmc_termination", this->builtin_node(value, false));
        else if (role_exists(type, BuiltinRole::SMCTermination))
          method("smc_termination", this->builtin_node(value, false));
        else
          find_builtin(type, BuiltinRole::MCMCTermination);
      }
      else if (target == "sequence")
      {
        // factors = [names]: the factors' numbers, as the engine reads them
        RecipeValue call = value;
        RecipeValue factor_numbers_value;
        if (kind_of(call) == "call")
        {
          const RecipeValue &arguments = call["arguments"];
          std::vector<std::string> keys = arguments.names();
          RecipeValue kept = RecipeValue::list();
          for (size_t k = 0; k < arguments.size(); ++k)
          {
            if (k < keys.size() && keys[k] == "factors")
              factor_numbers_value = this->factor_numbers_of(arguments[k]);
            else
              kept.push_back(arguments[k], k < keys.size() ? keys[k] : "");
          }
          call.set("arguments", kept);
        }
        RecipeValue node = this->builtin_node(call, true);
        if (!factor_numbers_value.is_null())
        {
          RecipeValue parameters = node["parameters"];
          parameters.push_back(factor_numbers_value, "factors");
          node.set("parameters", parameters);
        }
        // as the SMC sequence and as the EnK sequence (EnKI refuses the types
        // it cannot use when it reads it)
        method("smc_sequence", node);
        method("enk_sequence", node);
      }
      else if (target == "target")
        method("adaptive_target", this->builtin_node(value, false));
      else if (target == "resampling")
        method("adaptive_resampling", this->builtin_node(value, false));
      else if (target == "shifter")
        method("enk_shifter", this->builtin_node(value, false));
      else if (target == "weights")
      {
        // the probabilities of choosing each move: a function, or numbers
        arma::mat numbers;
        if (!this->is_function_value(value) && this->numeric(value, numbers))
        {
          std::string list;
          for (arma::uword k = 0; k < numbers.n_elem; ++k)
            list = list + (k ? ", " : "") + number_text(numbers(k));
          std::string function = "ilike_generated_weights_" + std::to_string(++this->generated_count);
          this->generated_code += "NumericVector " + function + "()\n{\n  return NumericVector::create(" + list + ");\n}\n\n";
          method("mcmc_weights", this->placeholder(function, "mcmc_weights"));
        }
        else
          method("mcmc_weights", this->function_placeholder(value, "mcmc_weights", "weights"));
      }
      else if (target == "proposal")
        is_proposals.push_back(this->is_proposal(value));
      else if (target == "transition_proposal")
        transition_proposals.push_back(this->transition_proposal(value));
      else
      {
        // a run-time setting (particles = 1000, …), read by the R function
        // that runs the algorithm (dev/ilike-language.md 4.7)
        if (settings.contains(target))
          this->error(statement, "\"" + target + "\" is set more than once.");
        settings.push_back(this->value_of(value), target);
      }
    }

    // the order of the moves
    std::vector<const Move*> ordered;
    if (has_sweep)
    {
      std::set<std::string> used;
      for (const std::string &n : sweep)
      {
        const Move *found = nullptr;
        for (const Move &m : moves)
          if (m.name == n) found = &m;
        if (found == nullptr)
          this->error(section, "sweep: there is no move called " + n + ".");
        ordered.push_back(found);
        used.insert(n);
      }
      for (const Move &m : moves)
        if (!used.count(m.name))
          this->error(section, "move " + m.name + " is not in the sweep.");
    }
    else
      for (const Move &m : moves)
        ordered.push_back(&m);

    RecipeValue algorithm = RecipeValue::list();
    std::vector<double> order;
    std::map<std::string, RecipeValue> by_kind;
    for (const Move *m : ordered)
    {
      order.push_back(m->code);
      std::string key = m->prefix + "_proposal";
      if (!by_kind.count(key))
        by_kind[key] = RecipeValue::list();
      by_kind[key].push_back(m->entry);
    }
    algorithm.push_back(order.empty() ? RecipeValue() : RecipeValue(order), "order_of_mcmc");
    for (const auto &kv : by_kind)
      algorithm.push_back(kv.second, kv.first);
    if (is_proposals.size() > 0)
      algorithm.push_back(is_proposals, "is_proposal");
    if (transition_proposals.size() > 0)
      algorithm.push_back(transition_proposals, "transition_proposal");
    if (methods.size() > 0)
      algorithm.push_back(methods, "method");
    algorithm.push_back(settings, "settings");
    this->algorithms.push_back(algorithm, name);
  }

  // transition_proposal = density(simulate = cpp::f, log = cpp::g), or a
  // built-in random walk on a state variable
  RecipeValue transition_proposal(const RecipeValue &value)
  {
    RecipeValue entry = RecipeValue::list();
    if (kind_of(value) == "call" && value["function"].as_string() == "density")
    {
      Arguments arguments = this->split_arguments(value);
      this->allow_only(value, arguments, {"simulate", "log"}, 0);
      if (!arguments.has("simulate") || !arguments.has("log"))
        this->error(value, "a transition proposal needs density(simulate = cpp::f, log = cpp::g).");
      entry.push_back(this->function_placeholder(arguments.named["simulate"], "simulate_transition_proposal", "simulate"),
                      "simulate_transition_proposal");
      entry.push_back(this->function_placeholder(arguments.named["log"], "evaluate_log_transition_proposal", "log"),
                      "evaluate_log_transition_proposal");
      return entry;
    }
    if (kind_of(value) != "call")
      this->error(value, "a transition proposal is density(simulate = cpp::f, log = cpp::g) or a built-in random walk.");
    std::string type = this->builtin_name(value);
    if (!role_exists(type, BuiltinRole::TransitionProposal))
      find_builtin(type, BuiltinRole::TransitionProposal);
    entry.push_back(this->builtin_node(value, true), "transition_proposal");
    return entry;
  }

  // proposal = norm(theta, 0, 3), or density(simulate = cpp::f, log = cpp::g)
  RecipeValue is_proposal(const RecipeValue &value)
  {
    RecipeValue entry = RecipeValue::list();
    if (kind_of(value) == "call" && value["function"].as_string() == "density")
    {
      Arguments arguments = this->split_arguments(value);
      this->allow_only(value, arguments, {"simulate", "log"}, 0);
      if (!arguments.has("simulate") || !arguments.has("log"))
        this->error(value, "an importance proposal written by the user needs density(simulate = cpp::f, log = cpp::g).");
      entry.push_back(this->function_placeholder(arguments.named["simulate"], "simulate_is_proposal", "simulate"),
                      "simulate_is_proposal");
      entry.push_back(this->function_placeholder(arguments.named["log"], "evaluate_log_is_proposal", "log"),
                      "evaluate_log_is_proposal");
      return entry;
    }
    if (kind_of(value) != "call" || !role_exists(this->builtin_name(value), BuiltinRole::ISProposal))
      this->error(value, "proposal is a built-in distribution (e.g. norm(theta, 0, 3)) or density(simulate = cpp::f, log = cpp::g).");
    entry.push_back(this->builtin_node(value, true), "is_proposal");
    return entry;
  }

};

}

RecipeValue lower_ilike_language(const RecipeValue &trees)
{
  Lowering lowering(trees);
  return lowering.run();
}

// The signature "type name(arguments)" of the definition of a C++ function
// in code (comments and strings skipped), or "" if there is none.
std::string find_cpp_function_signature(const std::string &code, const std::string &name)
{
  size_t n = code.size();
  size_t i = 0;
  size_t statement_start = 0;
  while (i < n)
  {
    char c = code[i];
    if (c == '/' && i + 1 < n && code[i + 1] == '/')
    {
      while (i < n && code[i] != '\n') ++i;
      statement_start = i;
      continue;
    }
    if (c == '/' && i + 1 < n && code[i + 1] == '*')
    {
      size_t end = code.find("*/", i + 2);
      i = (end == std::string::npos) ? n : end + 2;
      statement_start = i;
      continue;
    }
    if (c == '"' || c == '\'')
    {
      ++i;
      while (i < n && code[i] != c)
      {
        if (code[i] == '\\') ++i;
        ++i;
      }
      ++i;
      continue;
    }
    if (c == '#')
    {
      // a preprocessor line
      while (i < n && code[i] != '\n') ++i;
      statement_start = i;
      continue;
    }
    if (c == ';' || c == '{' || c == '}')
    {
      ++i;
      statement_start = i;
      continue;
    }
    if (std::isalpha((unsigned char)c) || c == '_')
    {
      size_t start = i;
      while (i < n && (std::isalnum((unsigned char)code[i]) || code[i] == '_')) ++i;
      bool preceded_by_member = start > 0 && (code[start - 1] == '.' || code[start - 1] == ':' ||
                                              (code[start - 1] == '>' && start > 1 && code[start - 2] == '-'));
      if (code.compare(start, i - start, name) == 0 && i - start == name.size() && !preceded_by_member)
      {
        size_t j = i;
        while (j < n && std::isspace((unsigned char)code[j])) ++j;
        if (j < n && code[j] == '(')
        {
          // the matching parenthesis, then "{" for a definition
          int depth = 0;
          size_t k = j;
          for (; k < n; ++k)
          {
            if (code[k] == '(') ++depth;
            else if (code[k] == ')' && --depth == 0) break;
          }
          size_t m = k + 1;
          while (m < n && std::isspace((unsigned char)code[m])) ++m;
          // something before the name (its return type) is needed
          std::string before = code.substr(statement_start, start - statement_start);
          bool has_type = before.find_first_not_of(" \t\r\n") != std::string::npos;
          if (k < n && m < n && code[m] == '{' && has_type)
          {
            std::string signature = code.substr(statement_start, k + 1 - statement_start);
            size_t first = signature.find_first_not_of(" \t\r\n");
            return signature.substr(first);
          }
        }
      }
      continue;
    }
    ++i;
  }
  return "";
}

}
