// Lowering the ilike language (dev/ilike-language.md; dev/plan.md 3.14, I3b):
// from the syntax trees of the files of a recipe to the recipe the engine
// reads — the same tree the old parser makes — with each user function
// marked by a placeholder {ilike_function, role} for the binding to compile
// and replace by the compiled function (role: the old block name, which fixes
// the function's signature). Small functions that the old format had users
// write (the name of a data variable, a list of factor numbers, literal data)
// are generated here as C++ code.

#include <functional>
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

  // values given for constants when compiling (compile(model_parameter_list)),
  // for what must be known then (the elements of moves on elements)
  void set_constant_values(const RecipeValue &values) { this->constant_values = values; }

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

    // "." in names: only in expressions in R (now that the imports' languages are known)
    for (size_t t = 0; t < this->trees.size(); ++t)
    {
      this->file = this->trees[t]["file"].as_string();
      this->check_dotted_names(this->trees[t]["sections"]);
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
      {
        this->finish_model_statements();
        this->finish_mfds();
      }
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
    if (this->mfds_entry.size() > 0)
      this->recipe.push_back(this->mfds_entry, "mfds");

    RecipeValue output = RecipeValue::list();
    output.push_back(this->recipe, "recipe");
    output.push_back(this->algorithms, "algorithms");
    output.push_back(this->constants, "constants");
    output.push_back(RecipeValue(this->declared_constants), "declared_constants");
    output.push_back(RecipeValue(this->user_code), "user_code");
    std::string all_generated = this->generated_code;
    if (!this->model_code.empty())
      all_generated += this->model_header + this->model_code + this->deferred_node_code + this->deferred_function_code;
    all_generated += this->mfds_code;
    output.push_back(RecipeValue(all_generated), "generated_code");
    output.push_back(RecipeValue(std::vector<std::string>(this->model_constants_used.begin(), this->model_constants_used.end())),
                     "model_constants");
    output.push_back(RecipeValue(std::vector<std::string>(this->element_constants_used.begin(), this->element_constants_used.end())),
                     "element_constants");
    output.push_back(this->foreign_code, "foreign_code");
    output.push_back(this->foreign_files, "foreign_files");
    output.push_back(RecipeValue(this->uses_foreign), "uses_foreign");
    output.push_back(RecipeValue(this->element_file), "element_file");
    return output;
  }

private:

  RecipeValue trees;
  std::string file;
  std::string element_file; // the file named in the generated code's messages
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
  RecipeValue constant_values;
  std::set<std::string> model_constants_used;
  std::set<std::string> element_constants_used; // (in moves on elements: values taken when compiling)
  // MFDS (dev/plan.md 3.16): the likelihood factors it can simulate, those
  // it cannot (with the reason), and the priors that are not Gaussian
  struct MfdsPart { std::string kind; RecipeValue call; std::string name; }; // kind: linear, nonlinear, sbi
  std::vector<MfdsPart> mfds_parts;
  std::vector<std::string> mfds_reasons;
  std::vector<std::string> mfds_prior_reasons;
  std::vector<std::string> mfds_log_variables; // simple lnorm/mvlnorm priors
  bool mfds_sbi_needs_direct = false;           // an sbi factor in another language, or with summaries
  RecipeValue mfds_entry = RecipeValue::list();
  std::string mfds_code;                         // emitted after the model's code
  // plates (I4b): the parameters numbered for the plates' lists of reads, and
  // the parameters each generated factor reads (by factor number)
  std::map<std::string, size_t> plate_variable_ids;
  std::vector<std::string> plate_variable_names;
  std::map<size_t, std::set<std::string>> generated_factor_reads;
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
    std::string op = statement.contains("operator") ? statement["operator"].as_string() : "";
    if (op == "+=" || op == "<-")
      this->error(statement, "\"" + op + "\" is only allowed in the model section, not in " + section + ".");
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

  // the call (or name) an expression is built on, through $ and [ ]: the
  // r::f(…) of r::f(…)$y[1]
  static const RecipeValue &foreign_base(const RecipeValue &expression)
  {
    const RecipeValue *current = &expression;
    while (kind_of(*current) == "member" || kind_of(*current) == "index")
      current = &(*current)["target"];
    return *current;
  }

  // TRUE for r::f, py::f, jl::f, alias::f (of an imported file), a call to
  // one, or $ / [ ] applied to such a call (r::read.csv("y.csv")$y)
  bool is_foreign(const RecipeValue &expression) const
  {
    const RecipeValue &base = foreign_base(expression);
    std::string kind = kind_of(base);
    if (kind != "name" && kind != "call")
      return false;
    if (kind == "name" && &base != &expression)
      return false; // (r::f$y: $ needs a call)
    std::string name_space = base["namespace"].as_string();
    return !name_space.empty() && name_space != "cpp" && name_space != "ilike" && this->foreign_languages.count(name_space);
  }

  // names may contain "." only in expressions in R (r::read.csv, or a call
  // through an alias of an imported .R file), as R allows; elsewhere (C++,
  // Python, Julia, ilike's own names) it is an error
  void check_dotted_names(const RecipeValue &node) const
  {
    if (!node.is_list())
      return;
    const RecipeValue &base = foreign_base(node);
    std::string kind = kind_of(base);
    if ((kind == "call" || kind == "name") && !base["namespace"].as_string().empty())
    {
      auto found = this->foreign_languages.find(base["namespace"].as_string());
      if (found != this->foreign_languages.end() && found->second == "r")
        return;
    }
    static const std::vector<std::string> fields = {"name", "function", "index", "alias", "keyword"};
    for (const std::string &field : fields)
      if (node.contains(field) && node[field].is_string() && node[field].as_string().find('.') != std::string::npos)
        this->error(node, "\"" + node[field].as_string() + "\": names cannot contain \".\" (except in calls to R, e.g. r::read.csv).");
    std::vector<std::string> keys = node.names();
    for (size_t k = 0; k < node.size(); ++k)
    {
      if (k < keys.size() && keys[k].find('.') != std::string::npos)
        this->error(node, "\"" + keys[k] + "\": names cannot contain \".\" (except in calls to R, e.g. r::read.csv).");
      this->check_dotted_names(node[k]);
    }
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
    {
      // parameters$theta, data$y: the objects of the role (the binding reads
      // these as its arguments); otherwise a member of a value in R, written
      // [["y"]] so that the binding does not take y for a variable
      const RecipeValue &target = expression["target"];
      if (kind_of(target) == "name" && target["namespace"].as_string().empty() &&
          special.count(target["name"].as_string()))
        return this->r_text(target) + "$" + expression["name"].as_string();
      return this->r_text(target) + "[[" + quoted(expression["name"].as_string()) + "]]";
    }
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
    std::string language = this->foreign_languages.at(foreign_base(expression)["namespace"].as_string());
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

  // a list of variables [a, b, c] as the first argument of a move's built-in
  // proposal: the engine's stacked name "a;b;c" (C203, C205), the column
  // vector of their elements in that order
  std::string stacked_variable(const RecipeValue &list, const std::string &type)
  {
    const RecipeValue &elements = list["elements"];
    if (elements.size() == 0)
      this->error(list, type + ": the list of variables is empty.");
    std::string stacked;
    std::set<std::string> seen;
    for (size_t k = 0; k < elements.size(); ++k)
    {
      if (kind_of(elements[k]) != "name" || this->is_function_reference(elements[k]))
        this->error(elements[k], type + ": each element of the list is the name of a variable.");
      std::string name = elements[k]["name"].as_string();
      if (!seen.insert(name).second)
        this->error(elements[k], type + ": " + name + " appears twice in the list of variables.");
      stacked = stacked + (stacked.empty() ? "" : ";") + name;
    }
    return stacked;
  }

  // a built-in call as the old recipe's node {type, variables, parameters}:
  // the variable is the given one, or the first positional argument (or the
  // argument "variable"); with_variable = false for settings without one
  RecipeValue builtin_node(const RecipeValue &call, bool with_variable, const std::string &variable = "",
                           bool lists = false)
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
        if (kind_of(arguments[i]) == "vector")
        {
          // several variables moved as one block (C203): the stacked name
          if (!lists)
            this->error(arguments[i], type + ": a list of variables can be moved together by the proposal of a move "
                        "(e.g. metropolis(mvnorm_rw([a, b], cov))), not here.");
          the_variable = this->stacked_variable(arguments[i], type);
          variable_taken = true;
          continue;
        }
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

  std::string generate_factor_index_numbers(const std::vector<size_t> &numbers)
  {
    std::string list;
    for (size_t n : numbers)
      list = list + (list.empty() ? "" : ", ") + std::to_string(n);
    std::string name = "ilike_generated_factor_index_" + std::to_string(++this->generated_count);
    this->generated_code += "NumericVector " + name + "()\n{\n  return NumericVector::create(" + list + ");\n}\n\n";
    return name;
  }

  // the factors a move on `variables` (a built-in proposal's, separated by ;
  // or ,) must evaluate: empty when that is every factor (or not known)
  std::vector<size_t> automatic_factor_index(const std::string &variables)
  {
    std::vector<size_t> numbers;
    if (variables.empty() || this->generated_factor_reads.empty())
      return numbers;
    std::set<std::string> moved;
    std::string current;
    for (char c : variables + ";")
    {
      if (c == ';' || c == ',')
      {
        if (!current.empty()) moved.insert(current);
        current.clear();
      }
      else if (c != ' ')
        current += c;
    }
    size_t total = this->factor_entries.size();
    for (size_t n = 1; n <= total; ++n)
    {
      auto found = this->generated_factor_reads.find(n);
      bool include = (found == this->generated_factor_reads.end());
      if (!include)
        for (const std::string &v : moved)
          if (found->second.count(v))
            include = true;
      if (include)
        numbers.push_back(n);
    }
    if (numbers.size() == total)
      numbers.clear();
    return numbers;
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
      if (keyword == "for" || op == "<-" || op == "+=" || (op == "~" && keyword.empty() && this->generate))
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
        std::string family = this->builtin_name(value);
        if (family == "lnorm" || family == "mvlnorm")
          this->mfds_log_variables.push_back(name);
        continue;
      }
      if (op != "=")
        this->error(statement, "a model statement is \"variable ~ distribution(…)\", \"prior name = …\" or \"factor name = …\".");

      if (keyword.empty() && name == "time")
      {
        // the time index of a state-space model: time = ssm(…)
        this->mfds_reasons.push_back("the model is a state-space model");
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
        this->mfds_prior_reasons.push_back("the prior " + name + " is written as a function");
        continue;
      }
      if (keyword == "transition")
      {
        this->mfds_reasons.push_back("the model is a state-space model");
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
        {
          this->add_factor(statement, name, this->sbi_factor(value));
          Arguments sbi_arguments = this->split_arguments(value);
          if (!this->is_function_reference(sbi_arguments.named["simulate"]) ||
              sbi_arguments.named["simulate"]["namespace"].as_string() != "cpp" || sbi_arguments.has("summaries"))
            this->mfds_sbi_needs_direct = true;
          this->mfds_parts.push_back({"sbi", value, name});
        }
        else if (function == "linear_gaussian")
        {
          this->add_factor(statement, name, this->gaussian_data_model(value, true));
          this->mfds_parts.push_back({"linear", value, name});
        }
        else if (function == "nonlinear_gaussian")
        {
          this->add_factor(statement, name, this->gaussian_data_model(value, false));
          this->mfds_parts.push_back({"nonlinear", value, name});
        }
        else if (!function.empty() && role_exists(function, BuiltinRole::AlgorithmicLikelihood))
        {
          this->add_factor(statement, name, this->algorithmic_likelihood(value));
          this->mfds_reasons.push_back("the factor " + name + " is estimated by a filter");
        }
        else if (!function.empty() && function != "density")
        {
          if (role_exists(function, BuiltinRole::LinearGaussianDataModel))
          {
            RecipeValue factor = RecipeValue::list();
            factor.push_back(this->builtin_node(value, true), "linear_gaussian_data_model");
            this->add_factor(statement, name, factor);
            this->mfds_reasons.push_back("the factor " + name + " (linear_gaussian_data_model) can be simulated "
                                         "when written linear_gaussian(…)");
          }
          else
            this->error(value, "unknown kind of factor " + function + " (factors: cpp::f, density(…), sbi(…), "
                        "linear_gaussian(…), nonlinear_gaussian(…), pf(…), kf(…), enkf(…)).");
        }
        else
        {
          this->add_factor(statement, name, this->density_functions(value, false));
          this->mfds_reasons.push_back("the factor " + name + " is a log density written as a function");
        }
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
          if (keyword == "for" || op == "<-" || op == "+=")
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
    if (op == "+=")
    {
      // a term of a named factor, not tied to a variable: factor name += expression
      if (keyword != "factor")
        this->error(statement, "\"+=\" adds a term to a named factor: factor name += expression.");
      if (kind_of(statement["target"]) != "name")
        this->error(statement, "a factor's name has no index (factor name += expression; the loop gives the terms).");
      ModelStatement st;
      st.op = op;
      st.node = statement;
      st.loops = loops;
      st.variable = statement["target"]["name"].as_string();
      st.indices = RecipeValue::list();
      st.value = statement["value"];
      st.observed = true; // (a likelihood factor: not part of the prior, and no target to read)
      if (this->data_names.count(st.variable) || this->constant_names.count(st.variable))
        this->error(statement, st.variable + " is " + (this->data_names.count(st.variable) ? "data" : "a constant") +
                    ", so it cannot name a factor of terms (choose another name).");
      if (!this->model_factor_position.count(st.variable))
      {
        this->define(statement, st.variable);
        this->factor_entries.push_back(RecipeValue::list());
        this->factor_numbers[st.variable] = this->factor_entries.size();
        this->model_factor_position[st.variable] = this->factor_entries.size() - 1;
      }
      this->model_statements.push_back(st);
      return;
    }
    if (!keyword.empty() || (op != "~" && op != "<-"))
      this->error(statement, "inside a loop, a statement is \"y[i] ~ distribution(…)\", \"mu[i] <- expression\" "
                  "or \"factor name += expression\".");
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
    bool node = false;                  // in a deterministic node's function (data is a pointer, nullptr in the prior)
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
  // a scalar constant's value when compiling: the value given to compile(),
  // or the constants section's default
  bool constant_number(const std::string &name, double &value) const
  {
    if (!this->constant_names.count(name))
      return false;
    const RecipeValue *found = NULL;
    if (this->constant_values.is_list() && this->constant_values.contains(name))
      found = &this->constant_values[name];
    else if (this->constants.contains(name))
      found = &this->constants[name];
    if (found == NULL || !found->is_numeric())
      return false;
    std::vector<double> numbers = found->as_doubles();
    if (numbers.size() != 1)
      return false;
    value = numbers[0];
    return true;
  }

  bool is_literal_data(const std::string &name) const
  {
    for (const auto &item : this->literal_data)
      if (item.first == name)
        return true;
    return false;
  }

  std::string literal_data_matrix(const std::string &name, const RecipeValue &node)
  {
    for (const auto &item : this->literal_data)
      if (item.first == name)
      {
        std::string symbol = "ilike_data_" + name;
        if (this->model_header.find("static arma::mat " + symbol + "_make") == std::string::npos)
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
      // a node is also called by the prior, with no data: literal data are
      // then read from their values compiled in (data given by a function
      // cannot reach the prior, deterministic_call)
      if (context.node && this->is_literal_data(name))
        return "(data != nullptr ? (*data)[" + quoted(name) + "] : " + this->literal_data_matrix(name, node) + ")";
      return context.data + "[" + quoted(name) + "]";
    }
    if (this->is_deterministic(name))
      this->error(node, name + " is a deterministic node: use it with its indices.");
    return "parameters[" + quoted(name) + "]";
  }

  // "x[...] (file:line:column)", for messages about an element. The file is
  // a marker, replaced when the message is made by the file name R sets after
  // loading the model (ilike_model_set_file): the generated code, and so the
  // compile cache's key, then do not depend on where the file is (C208).
  std::string element_description(const std::string &name, const RecipeValue &node)
  {
    this->element_file = this->file;
    std::string where = "@ILIKE_FILE@";
    if (node.is_list() && node.contains("line"))
      where = where + ":" + std::to_string(node["line"].as_int()) + ":" + std::to_string(node["column"].as_int());
    return name + "[...] (" + where + ")";
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
      bool stochastic = false;
      for (size_t k = 0; k < indices.size(); ++k)
        if (!this->parameters_read(indices[k], context).empty())
          stochastic = true;
      std::string checked = std::string(stochastic ? "true" : "false") + ", " + quoted(this->element_description(name, expression));
      if (indices.size() == 1)
        return "ilike_model_at(" + base + ", " + this->index_code(indices[0], context) + ", " + checked + ")";
      if (indices.size() == 2)
        return "ilike_model_at(" + base + ", " + this->index_code(indices[0], context) + ", " +
               this->index_code(indices[1], context) + ", " + checked + ")";
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
        // (both double: a number such as 2 is written as an int literal)
        return "std::" + f + "((double)" + argument(0) + ", (double)" + argument(1) + ")";
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
      // (the nested list's type named: arma::mat({{0}, {0}}) is ambiguous for
      // g++, which can read the inner braces as braced scalars of a flat list)
      std::string code = "arma::mat(std::initializer_list<std::initializer_list<double>>{";
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
          // (the node's own indices are not parameters)
          const ModelStatement &node = this->deterministic_statement(name);
          std::set<std::string> inner;
          this->read_parameters(node.value, inner, visited);
          for (size_t k = 0; k < node.indices.size(); ++k)
            if (kind_of(node.indices[k]) == "name")
              inner.erase(node.indices[k]["name"].as_string());
          found.insert(inner.begin(), inner.end());
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
    // (inside a node's function, data is already a pointer)
    std::string data = context.prior ? std::string("nullptr") : (context.node ? std::string("data") : std::string("&data"));
    std::string call = "ilike_model_node_" + name + "(parameters, " + data;
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

  // a statement of a likelihood factor: an observed variable's "~", or a
  // "factor name += expression" term
  bool likelihood_term(const ModelStatement &st) const
  {
    return (st.op == "~" && st.observed) || st.op == "+=";
  }

  // one term's log value: the log density of a "~", or the expression of a "+="
  std::string term_code(const ModelStatement &st, const CodeContext &context)
  {
    if (st.op == "+=")
      return "(double)(" + this->scalar_code(st.value, context) + ")";
    return this->density_code(st, this->target_code(st, context, false), context);
  }

  std::string random_code(const ModelStatement &st, const CodeContext &context)
  {
    const RecipeValue &call = st.value;
    std::string f = call["function"].as_string();
    const RecipeValue &arguments = call["arguments"];
    static const std::map<std::string, std::string> samplers = {
      {"norm", "rnorm"}, {"lnorm", "rlnorm"}, {"gamma", "rgamma"}, {"unif", "runif"}, {"exp", "rexp"},
      {"pois", "(double)rpois"}, {"tnorm", "rtnorm"}, {"beta", "ilike_model_rbeta"},
      {"bern", "ilike_model_rbern"}, {"binom", "ilike_model_rbinom"}};
    if (f == "mvnorm")
      return "arma::mat(rmvnorm(rng, arma::colvec(arma::vectorise(" + this->matrix_code(arguments[0], context) + ")), " +
             this->matrix_code(arguments[1], context) + "))";
    auto found = samplers.find(f);
    if (found == samplers.end())
      this->error(call, "the parameter " + st.variable + " has the distribution " + f +
                  ", which cannot be simulated in a prior.");
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
    if (st.indices.size() > 2)
      this->error(st.node, "at most two indices.");
    if (writable)
    {
      if (st.indices.size() == 1)
        return "(" + base + ")(" + this->index_code(st.indices[0], context) + ")";
      return "(" + base + ")(" + this->index_code(st.indices[0], context) + ", " + this->index_code(st.indices[1], context) + ")";
    }
    // a target read (in a log density): checked, as other elements
    bool stochastic = false;
    for (size_t k = 0; k < st.indices.size(); ++k)
      if (!this->parameters_read(st.indices[k], context).empty())
        stochastic = true;
    std::string checked = std::string(stochastic ? "true" : "false") + ", " + quoted(this->element_description(st.variable, st.node));
    if (st.indices.size() == 1)
      return "ilike_model_at(" + base + ", " + this->index_code(st.indices[0], context) + ", " + checked + ")";
    return "ilike_model_at(" + base + ", " + this->index_code(st.indices[0], context) + ", " +
           this->index_code(st.indices[1], context) + ", " + checked + ")";
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

  // --- gradients (I4e): code generated from each term ---------------------------
  //
  // Each term passes an adjoint (d term / d subexpression, a C++ expression)
  // down its expression; at a parameter's element the adjoint is added to
  // that variable's gradient, if it is one of the variables asked for
  // (ilike_G, a set of variables and their gradients: several variables are
  // computed in one pass over the terms, C204). Deterministic
  // nodes have gradient functions taking the adjoint. Subexpressions that read
  // no parameters are skipped. dev/plan.md, "Gradients of model statements".

  size_t gradient_temporaries = 0;

  static std::string indented_code(const std::string &code)
  {
    std::string out, line;
    std::stringstream lines(code);
    while (std::getline(lines, line))
      out += "  " + line + "\n";
    return out;
  }

  std::string node_gradient_call(const std::string &name, const RecipeValue &indices, const std::string &adjoint,
                                 const CodeContext &context)
  {
    std::string data = context.prior ? std::string("nullptr") : (context.node ? std::string("data") : std::string("&data"));
    std::string call = "ilike_model_node_gradient_" + name + "(ilike_G, " + adjoint + ", parameters, " + data;
    for (size_t k = 0; k < indices.size(); ++k)
      call += ", (int)(" + this->scalar_code(indices[k], context) + ")";
    return call + ");\n";
  }

  std::string gradient_code(const RecipeValue &expression, const std::string &adjoint, const CodeContext &context)
  {
    if (this->parameters_read(expression, context).empty())
      return "";
    std::string kind = kind_of(expression);
    if (kind == "name")
    {
      const std::string name = expression["name"].as_string();
      if (this->is_deterministic(name))
        return this->node_gradient_call(name, RecipeValue::list(), adjoint, context);
      return "if (arma::mat *ilike_gv = ilike_G.find(" + quoted(name) + ")) (*ilike_gv)(0) += " + adjoint + ";\n";
    }
    if (kind == "index")
    {
      const std::string name = expression["target"]["name"].as_string();
      const RecipeValue &indices = expression["indices"];
      if (this->is_deterministic(name))
        return this->node_gradient_call(name, indices, adjoint, context);
      if (this->constant_names.count(name) || this->data_names.count(name))
        return ""; // (an index that reads parameters is piecewise constant in them: no contribution)
      bool stochastic = false;
      for (size_t k = 0; k < indices.size(); ++k)
        if (!this->parameters_read(indices[k], context).empty())
          stochastic = true;
      std::string where = std::string(stochastic ? "true" : "false") + ", " + quoted(this->element_description(name, expression));
      std::string element = this->index_code(indices[0], context);
      if (indices.size() == 2)
        element += ", " + this->index_code(indices[1], context);
      return "if (arma::mat *ilike_gv = ilike_G.find(" + quoted(name) + ")) ilike_model_add(*ilike_gv, " + element + ", " + adjoint + ", " + where + ");\n";
    }
    std::string t = "ilike_a" + std::to_string(++this->gradient_temporaries);
    std::string open = "{\n  const double " + t + " = " + adjoint + ";\n";
    if (kind == "unary")
      return open + indented_code(this->gradient_code(expression["operand"], "-" + t, context)) + "}\n";
    if (kind == "binary")
    {
      std::string op = expression["operator"].as_string();
      const RecipeValue &left = expression["left"], &right = expression["right"];
      std::string a = this->scalar_code(left, context), b = this->scalar_code(right, context);
      std::string inner;
      if (op == "+")
        inner = this->gradient_code(left, t, context) + this->gradient_code(right, t, context);
      else if (op == "-")
        inner = this->gradient_code(left, t, context) + this->gradient_code(right, "-" + t, context);
      else if (op == "*")
        inner = this->gradient_code(left, t + "*" + b, context) + this->gradient_code(right, t + "*" + a, context);
      else if (op == "/")
        inner = this->gradient_code(left, t + "/" + b, context) +
                this->gradient_code(right, "-" + t + "*" + a + "/(" + b + "*" + b + ")", context);
      else if (op == "^")
        inner = this->gradient_code(left, t + "*" + b + "*std::pow(" + a + ", " + b + " - 1.0)", context) +
                this->gradient_code(right, t + "*std::log(" + a + ")*std::pow(" + a + ", " + b + ")", context);
      else
        this->error(expression, "the operator " + op + " cannot be differentiated.");
      return open + indented_code(inner) + "}\n";
    }
    if (kind == "call")
    {
      std::string f = expression["function"].as_string();
      const RecipeValue &arguments = expression["arguments"];
      std::string u = this->scalar_code(arguments[0], context);
      static const std::map<std::string, std::string> derivative = {
        {"exp", "std::exp(U)"}, {"log", "(1.0/(U))"}, {"sqrt", "(0.5/std::sqrt(U))"},
        {"abs", "(double)(((U) > 0.0) - ((U) < 0.0))"}, {"lgamma", "R::digamma(U)"}, {"log1p", "(1.0/(1.0 + (U)))"},
        {"expm1", "std::exp(U)"}, {"sin", "std::cos(U)"}, {"cos", "(-std::sin(U))"},
        {"tanh", "(1.0 - std::pow(std::tanh(U), 2))"},
        {"inv_logit", "(1.0/(1.0 + std::exp(-(U))))*(1.0 - 1.0/(1.0 + std::exp(-(U))))"},
        {"ilogit", "(1.0/(1.0 + std::exp(-(U))))*(1.0 - 1.0/(1.0 + std::exp(-(U))))"},
        {"logit", "(1.0/((U)*(1.0 - (U))))"}};
      auto found = derivative.find(f);
      if (found != derivative.end())
      {
        std::string d = found->second;
        for (size_t at = d.find('U'); at != std::string::npos; at = d.find('U', at + u.size()))
          d.replace(at, 1, u);
        return open + indented_code(this->gradient_code(arguments[0], t + "*" + d, context)) + "}\n";
      }
      std::string v = this->scalar_code(arguments[1], context);
      if (f == "pow")
        return open + indented_code(this->gradient_code(arguments[0], t + "*" + v + "*std::pow(" + u + ", " + v + " - 1.0)", context) +
                                    this->gradient_code(arguments[1], t + "*std::log(" + u + ")*std::pow(" + u + ", " + v + ")", context)) + "}\n";
      if (f == "min" || f == "max")
      {
        // (the derivative of the argument chosen)
        std::string first = (f == "min") ? "((double)" + u + " <= (double)" + v + ")" : "((double)" + u + " >= (double)" + v + ")";
        return open + "  if " + first + "\n  {\n" + indented_code(indented_code(this->gradient_code(arguments[0], t, context))) +
               "  }\n  else\n  {\n" + indented_code(indented_code(this->gradient_code(arguments[1], t, context))) + "  }\n}\n";
      }
      this->error(expression, "the function " + f + " cannot be differentiated.");
    }
    this->error(expression, "this expression cannot be differentiated.");
  }

  // the gradient of one "~" term: the log density's partial derivatives with
  // respect to the target and the arguments, then the arguments' chain rule
  std::string density_gradient_code(const ModelStatement &st, const CodeContext &context)
  {
    const RecipeValue &call = st.value;
    std::string f = call["function"].as_string();
    const RecipeValue &arguments = call["arguments"];
    if (f == "mvnorm")
    {
      std::string code = "{\n  arma::colvec ilike_x = arma::vectorise(" + this->target_code(st, context, false) + ");\n"
                         "  arma::colvec ilike_m = arma::vectorise(" + this->matrix_code(arguments[0], context) + ");\n"
                         "  arma::mat ilike_S = " + this->matrix_code(arguments[1], context) + ";\n"
                         "  arma::colvec ilike_r = arma::solve(ilike_S, ilike_x - ilike_m, arma::solve_opts::likely_sympd);\n";
      auto is_parameter = [&](const RecipeValue &e) {
        return kind_of(e) == "name" && !this->data_names.count(e["name"].as_string()) &&
               !this->constant_names.count(e["name"].as_string()); };
      if (!st.observed)
        code += "  if (arma::mat *ilike_gv = ilike_G.find(" + quoted(st.variable) + ")) *ilike_gv -= arma::reshape(ilike_r, ilike_gv->n_rows, ilike_gv->n_cols);\n";
      if (is_parameter(arguments[0]))
        code += "  if (arma::mat *ilike_gv = ilike_G.find(" + quoted(arguments[0]["name"].as_string()) +
                ")) *ilike_gv += arma::reshape(ilike_r, ilike_gv->n_rows, ilike_gv->n_cols);\n";
      if (is_parameter(arguments[1]))
        code += "  if (ilike_G.find(" + quoted(arguments[1]["name"].as_string()) + ")) Rcpp::stop(" +
                quoted("The gradient with respect to the covariance of mvnorm (" + arguments[1]["name"].as_string() +
                       ") is not available yet.") + ");\n";
      return code + "}\n";
    }
    std::string x = "ilike_x", v[4] = {"ilike_v0", "ilike_v1", "ilike_v2", "ilike_v3"};
    std::string dx;
    std::vector<std::string> d(arguments.size());
    std::string prelude;
    if (f == "norm")
    {
      dx = "(-(" + x + " - " + v[0] + ")/(" + v[1] + "*" + v[1] + "))";
      d[0] = "((" + x + " - " + v[0] + ")/(" + v[1] + "*" + v[1] + "))";
      d[1] = "(-1.0/" + v[1] + " + (" + x + " - " + v[0] + ")*(" + x + " - " + v[0] + ")/(" + v[1] + "*" + v[1] + "*" + v[1] + "))";
    }
    else if (f == "lnorm")
    {
      std::string l = "(std::log(" + x + ") - " + v[0] + ")";
      dx = "(-1.0/" + x + " - " + l + "/(" + v[1] + "*" + v[1] + "*" + x + "))";
      d[0] = "(" + l + "/(" + v[1] + "*" + v[1] + "))";
      d[1] = "(-1.0/" + v[1] + " + " + l + "*" + l + "/(" + v[1] + "*" + v[1] + "*" + v[1] + "))";
    }
    else if (f == "gamma")
    {
      dx = "((" + v[0] + " - 1.0)/" + x + " - " + v[1] + ")";
      d[0] = "(std::log(" + v[1] + ") - R::digamma(" + v[0] + ") + std::log(" + x + "))";
      d[1] = "(" + v[0] + "/" + v[1] + " - " + x + ")";
    }
    else if (f == "exp")
    {
      dx = "(-" + v[0] + ")";
      d[0] = "(1.0/" + v[0] + " - " + x + ")";
    }
    else if (f == "unif")
    {
      d[0] = "(1.0/(" + v[1] + " - " + v[0] + "))";
      d[1] = "(-1.0/(" + v[1] + " - " + v[0] + "))";
    }
    else if (f == "tnorm")
    {
      prelude = "  const arma::colvec ilike_d = ilike_model_dtnorm_partials(" + x + ", " + v[0] + ", " + v[1] + ", " + v[2] + ", " + v[3] + ");\n";
      dx = "ilike_d[0]";
      for (size_t k = 0; k < 4; ++k)
        d[k] = "ilike_d[" + std::to_string(k + 1) + "]";
    }
    else if (f == "beta")
    {
      dx = "((" + v[0] + " - 1.0)/" + x + " - (" + v[1] + " - 1.0)/(1.0 - " + x + "))";
      d[0] = "(std::log(" + x + ") - R::digamma(" + v[0] + ") + R::digamma(" + v[0] + " + " + v[1] + "))";
      d[1] = "(std::log(1.0 - " + x + ") - R::digamma(" + v[1] + ") + R::digamma(" + v[0] + " + " + v[1] + "))";
    }
    else if (f == "pois")
      d[0] = "(" + x + "/" + v[0] + " - 1.0)"; // (the target is discrete)
    else if (f == "bern")
      d[0] = "(" + x + "/" + v[0] + " - (1.0 - " + x + ")/(1.0 - " + v[0] + "))";
    else if (f == "binom")
      d[1] = "(" + x + "/" + v[1] + " - (" + v[0] + " - " + x + ")/(1.0 - " + v[1] + "))"; // (n is discrete)
    else
      this->error(call, "the distribution " + f + " cannot be differentiated.");

    std::string code = "{\n  const double " + x + " = " + this->target_code(st, context, false) + ";\n";
    for (size_t k = 0; k < arguments.size() && k < 4; ++k)
      code += "  const double " + v[k] + " = " + this->scalar_code(arguments[k], context) + ";\n";
    code += prelude;
    if (!st.observed && !dx.empty())
    {
      // the target: an element of a parameter (its indices are loop indices or numbers)
      std::string element = "0";
      if (st.indices.size() >= 1)
        element = this->index_code(st.indices[0], context);
      if (st.indices.size() == 2)
        element += ", " + this->index_code(st.indices[1], context);
      code += "  if (arma::mat *ilike_gv = ilike_G.find(" + quoted(st.variable) + ")) ilike_model_add(*ilike_gv, " + element + ", " + dx + ", false, " +
              quoted(this->element_description(st.variable, st.node)) + ");\n";
    }
    for (size_t k = 0; k < arguments.size(); ++k)
      if (!d[k].empty())
        code += indented_code(this->gradient_code(arguments[k], d[k], context));
    return code + "}\n";
  }

  std::string term_gradient_code(const ModelStatement &st, const CodeContext &context)
  {
    if (st.op == "+=")
      return this->gradient_code(st.value, "1.0", context);
    return this->density_gradient_code(st, context);
  }

  // the gradient functions of a factor made of these statements (one term per
  // statement and loop index): f_many(variables, gradients, parameters[,
  // data]), the gradients of several variables in one pass over the terms
  // (C204), and f(variable, parameters[, data]) for one variable; for the
  // prior also forms with data (for its plate). "" if they cannot be
  // generated (the factor then has no gradient, as before).
  std::string generate_gradient(const std::string &suffix, const std::vector<const ModelStatement*> &statements, bool prior)
  {
    std::string name = "ilike_model_gradient_" + suffix;
    try
    {
      std::string body = "ilike_model_gradients ilike_G(variables, gradients, parameters);\ntry\n{\n";
      for (const ModelStatement* st : statements)
      {
        std::set<std::string> reads = this->statements_read(std::vector<const ModelStatement*>(1, st));
        if (reads.empty())
          continue;
        std::string guard;
        for (const std::string &r : reads)
          guard += (guard.empty() ? "" : " || ") + std::string("ilike_G.find(") + quoted(r) + ")";
        CodeContext context;
        context.prior = prior;
        for (const ModelLoop &loop : st->loops) context.indices.insert(loop.index);
        std::string term = this->term_gradient_code(*st, context);
        CodeContext fresh;
        fresh.prior = prior;
        body += "  if (" + guard + ")\n  {\n" + indented_code(indented_code(this->loops_code(*st, fresh, term, ""))) + "  }\n";
      }
      body += "}\ncatch (const ilike_model_outside &)\n{\n  for (arma::mat &g : gradients) g.fill(arma::datum::nan); // (density 0 at this point)\n}\n";
      std::string many = "(const std::vector<std::string> &variables, std::vector<arma::mat> &gradients, const Parameters &parameters" +
                         std::string(prior ? ")" : ", const Data &data)");
      std::string code = "void " + name + "_many" + many + "\n{\n" + indented_code(body) + "}\n\n";
      std::string signature = prior ? "(const std::string &variable, const Parameters &parameters)"
                                    : "(const std::string &variable, const Parameters &parameters, const Data &data)";
      code += "arma::mat " + name + signature + "\n{\n  std::vector<std::string> variables(1, variable);\n  std::vector<arma::mat> gradients;\n  " +
              name + "_many(variables, gradients, parameters" + (prior ? "" : ", data") + ");\n  return gradients[0];\n}\n\n";
      if (prior)
        code += "arma::mat " + name + "_plate(const std::string &variable, const Parameters &parameters, const Data &data)\n{\n"
                "  return " + name + "(variable, parameters);\n}\n\n"
                "void " + name + "_many_plate(const std::vector<std::string> &variables, std::vector<arma::mat> &gradients, "
                "const Parameters &parameters, const Data &data)\n{\n  " + name + "_many(variables, gradients, parameters);\n}\n\n";
      this->deferred_function_code += code;
      return name;
    }
    catch (std::exception &e)
    {
      return "";
    }
  }

  // --- plates (I4b): what each term reads, and one term at a time ---------------

  size_t plate_variable_id(const std::string &name)
  {
    auto found = this->plate_variable_ids.find(name);
    if (found != this->plate_variable_ids.end())
      return found->second;
    size_t id = this->plate_variable_names.size();
    this->plate_variable_ids[name] = id;
    this->plate_variable_names.push_back(name);
    return id;
  }

  std::string read_code(const std::string &name, const std::string &row, const std::string &col)
  {
    return "ilike_model_read(reads, term, " + std::to_string(this->plate_variable_id(name)) + ", " + row + ", " + col + ");\n";
  }

  // the parameters an expression reads, other than the loop indices in scope
  std::set<std::string> parameters_read(const RecipeValue &expression, const CodeContext &context) const
  {
    std::set<std::string> found;
    this->read_parameters(expression, found);
    for (const std::string &index : context.indices)
      found.erase(index);
    return found;
  }

  // C++ recording what an expression reads, for the plate's structure
  // function: an element (row, column, 0-based) or a whole variable (-1, -1).
  // It mirrors scalar_code (and matrix_code when whole). An index that reads
  // parameters (stochastic indexing) cannot be evaluated without them, so the
  // indexed variable is then read whole.
  std::string reads_code(const RecipeValue &expression, const CodeContext &context, bool whole)
  {
    std::string kind = kind_of(expression);
    if (kind == "name")
    {
      if (!expression["namespace"].as_string().empty())
        return "";
      std::string name = expression["name"].as_string();
      if (context.indices.count(name) || this->constant_names.count(name) || this->data_names.count(name))
        return "";
      if (this->is_deterministic(name))
        return this->node_reads_call(name, RecipeValue::list(), context);
      return whole ? this->read_code(name, "-1", "-1") : this->read_code(name, "0", "0");
    }
    if (kind == "index")
    {
      const RecipeValue &target = expression["target"];
      std::string name = target["name"].as_string();
      const RecipeValue &indices = expression["indices"];
      std::string code;
      bool stochastic = false;
      for (size_t k = 0; k < indices.size(); ++k)
      {
        code += this->reads_code(indices[k], context, false);
        if (!this->parameters_read(indices[k], context).empty())
          stochastic = true;
      }
      if (this->is_deterministic(name))
      {
        if (stochastic)
          return code + this->node_reads_whole(name);
        return code + this->node_reads_call(name, indices, context);
      }
      if (this->constant_names.count(name) || this->data_names.count(name))
        return code;
      if (stochastic)
        return code + this->read_code(name, "-1", "-1");
      std::string row = this->index_code(indices[0], context);
      std::string col = (indices.size() == 2) ? this->index_code(indices[1], context) : std::string("0");
      return code + this->read_code(name, row, col);
    }
    if (!expression.is_list())
      return "";
    std::string code;
    for (size_t k = 0; k < expression.size(); ++k)
      if (expression[k].is_list())
        code += this->reads_code(expression[k], context, false);
    return code;
  }

  // every parameter a deterministic node reads, whole
  std::string node_reads_whole(const std::string &name)
  {
    std::set<std::string> found;
    this->read_parameters(this->deterministic_statement(name).value, found);
    const ModelStatement &st = this->deterministic_statement(name);
    for (size_t k = 0; k < st.indices.size(); ++k)
      found.erase(st.indices[k]["name"].as_string());
    std::string code;
    for (const std::string &p : found)
      code += this->read_code(p, "-1", "-1");
    return code;
  }

  std::string node_reads_call(const std::string &name, const RecipeValue &indices, const CodeContext &context)
  {
    const ModelStatement &st = this->deterministic_statement(name);
    if (indices.size() != st.indices.size())
      return this->node_reads_whole(name);
    std::string data = context.prior ? std::string("nullptr") : (context.node ? std::string("data") : std::string("&data"));
    std::string call = "ilike_model_node_reads_" + name + "(reads, term, " + data;
    for (size_t k = 0; k < indices.size(); ++k)
      call += ", (int)(" + this->scalar_code(indices[k], context) + ")";
    return call + ");\n";
  }

  // the reads of a statement: its target (a parameter) and its distribution's arguments
  std::string statement_reads_code(const ModelStatement &st, const CodeContext &context)
  {
    if (st.op == "+=")
      return this->reads_code(st.value, context, false);
    std::string code;
    bool mvnorm = st.value["function"].as_string() == "mvnorm";
    if (!st.observed)
    {
      if (mvnorm || st.indices.size() == 0)
        code += mvnorm ? this->read_code(st.variable, "-1", "-1") : this->read_code(st.variable, "0", "0");
      else
        code += this->read_code(st.variable, this->index_code(st.indices[0], context),
                                (st.indices.size() == 2) ? this->index_code(st.indices[1], context) : std::string("0"));
    }
    const RecipeValue &arguments = st.value["arguments"];
    for (size_t k = 0; k < arguments.size(); ++k)
      code += this->reads_code(arguments[k], context, mvnorm);
    return code;
  }

  // a plate needs a number of terms that does not depend on the parameters
  bool plate_possible(const std::vector<const ModelStatement*> &statements) const
  {
    for (const ModelStatement* st : statements)
      for (const ModelLoop &loop : st->loops)
      {
        std::set<std::string> found;
        this->read_parameters(loop.from, found);
        this->read_parameters(loop.to, found);
        for (const ModelLoop &other : st->loops) found.erase(other.index);
        if (!found.empty())
          return false;
      }
    return true;
  }

  // the structure and term functions of the plate of a factor whose terms are
  // the statements' loops, in this order (the order of the factor's sum)
  void generate_plate(const std::string &suffix, const std::vector<const ModelStatement*> &statements, bool prior)
  {
    size_t width = 1;
    for (const ModelStatement* st : statements)
      width = std::max(width, st->loops.size() + 1);
    std::string structure = "Data ilike_model_plate_" + suffix + "(const Data &data)\n{\n"
                            "  std::vector<double> terms, reads;\n  double term = 0.0;\n"
                            "  Parameters parameters; // (no index read here depends on the parameters)\n";
    std::string term = "double ilike_model_term_" + suffix +
                       "(const Parameters &parameters, const Data &data, const arma::mat &terms, size_t term)\n{\n"
                       "  try\n  {\n"
                       "  switch ((int)terms(term, 0))\n  {\n";
    for (size_t s = 0; s < statements.size(); ++s)
    {
      const ModelStatement &st = *statements[s];
      CodeContext context;
      context.prior = prior;
      for (const ModelLoop &loop : st.loops) context.indices.insert(loop.index);
      std::string body = "terms.push_back(" + std::to_string(s) + ");\n";
      for (const ModelLoop &loop : st.loops) body += "terms.push_back((double)" + loop.index + ");\n";
      for (size_t k = st.loops.size() + 1; k < width; ++k) body += "terms.push_back(0.0);\n";
      body += this->statement_reads_code(st, context);
      body += "term += 1.0;\n";
      CodeContext fresh;
      fresh.prior = prior;
      structure += this->loops_code(st, fresh, body, "  ");

      term += "    case " + std::to_string(s) + ":\n    {\n";
      for (size_t k = 0; k < st.loops.size(); ++k)
        term += "      int " + st.loops[k].index + " = (int)terms(term, " + std::to_string(k + 1) + ");\n";
      term += "      return " + this->term_code(st, context) + ";\n    }\n";
    }
    structure += "  Data result;\n  result[\"terms\"] = ilike_model_table(terms, " + std::to_string(width) + ");\n"
                 "  result[\"reads\"] = ilike_model_table(reads, 4);\n  return result;\n}\n\n";
    term += "  }\n  }\n" + std::string("  catch (const ilike_model_outside &)\n  {\n    return -std::numeric_limits<double>::infinity(); // (an index that reads parameters is outside its variable)\n  }\n") + "  return 0.0;\n}\n\n";
    this->deferred_function_code += structure + term;
  }

  // the parameters a factor's statements read (with the prior's targets)
  std::set<std::string> statements_read(const std::vector<const ModelStatement*> &statements) const
  {
    std::set<std::string> found;
    for (const ModelStatement* st : statements)
    {
      std::set<std::string> one;
      this->read_parameters(st->value, one);
      for (size_t i = 0; i < st->indices.size(); ++i) this->read_parameters(st->indices[i], one);
      for (const ModelLoop &loop : st->loops) { this->read_parameters(loop.from, one); this->read_parameters(loop.to, one); }
      for (const ModelLoop &loop : st->loops) one.erase(loop.index);
      if (!st->observed) one.insert(st->variable);
      found.insert(one.begin(), one.end());
    }
    return found;
  }

  RecipeValue plate_entry(const std::string &suffix, const std::string &factor_name)
  {
    RecipeValue plate = RecipeValue::list();
    plate.push_back(this->placeholder("ilike_model_plate_" + suffix, "plate_structure"), "structure");
    plate.push_back(this->placeholder("ilike_model_term_" + suffix, "plate_term"), "term");
    plate.push_back(RecipeValue(factor_name), "name");
    return plate;
  }

  // --- MFDS (dev/plan.md 3.16): a joint simulator of the likelihood, the
  // Gaussian prior's moments, and the exp of log-normal variables -------------

  // the data variables an expression reads, through the deterministic nodes it uses
  void data_read_by(const RecipeValue &node, std::set<std::string> &found, std::set<std::string> &nodes_seen) const
  {
    if (!node.is_list())
      return;
    if (kind_of(node) == "name" && node["namespace"].as_string().empty())
    {
      std::string name = node["name"].as_string();
      if (this->data_names.count(name))
        found.insert(name);
      if (this->is_deterministic(name) && !nodes_seen.count(name))
      {
        nodes_seen.insert(name);
        this->data_read_by(this->deterministic_statement(name).value, found, nodes_seen);
      }
    }
    for (size_t k = 0; k < node.size(); ++k)
      this->data_read_by(node[k], found, nodes_seen);
  }

  // a matrix of a Gaussian data model as C++: numbers, or a C++ function
  // (with or without the parameters); "" if it cannot be (another language)
  std::string mfds_matrix_code(const RecipeValue &expression)
  {
    arma::mat value;
    if (!this->is_function_value(expression) && this->numeric(expression, value))
    {
      std::string list;
      for (arma::uword k = 0; k < value.n_elem; ++k)
        list += (k ? ", " : "") + number_text(value(k));
      return "arma::mat(std::vector<double>{" + list + "}.data(), " + std::to_string(value.n_rows) + ", " +
             std::to_string(value.n_cols) + ")";
    }
    if (!this->is_function_reference(expression) || expression["namespace"].as_string() != "cpp")
      return "";
    std::string name = this->cpp_function(expression, "a matrix");
    return name + (this->functions_take_parameters({expression}) ? "(parameters)" : "()");
  }

  // the names of a data argument (y, or [y1, y2])
  static std::vector<std::string> mfds_names(const RecipeValue &expression)
  {
    std::vector<std::string> names;
    if (kind_of(expression) == "name")
      names.push_back(expression["name"].as_string());
    else if (kind_of(expression) == "vector")
      for (size_t k = 0; k < expression["elements"].size(); ++k)
        if (kind_of(expression["elements"][k]) == "name")
          names.push_back(expression["elements"][k]["name"].as_string());
    return names;
  }

  static std::string joined(const std::vector<std::string> &items, const std::string &separator)
  {
    std::string out;
    for (size_t k = 0; k < items.size(); ++k)
      out += (k ? separator : "") + items[k];
    return out;
  }

  void finish_mfds()
  {
    std::vector<std::string> reasons = this->mfds_reasons;
    std::vector<std::string> observed;
    std::vector<std::string> potentials;
    for (const ModelStatement &st : this->model_statements)
    {
      if (st.op == "~" && st.observed && std::find(observed.begin(), observed.end(), st.variable) == observed.end())
        observed.push_back(st.variable);
      if (st.op == "+=" && std::find(potentials.begin(), potentials.end(), st.variable) == potentials.end())
      {
        potentials.push_back(st.variable);
        reasons.push_back("the factor " + st.variable + " is made of += terms, which cannot be simulated");
      }
    }

    // the likelihood's simulator: generated, unless the likelihood is one
    // sbi factor, which MFDS uses directly (in any language, with summaries)
    size_t likelihoods = observed.size() + this->mfds_parts.size();
    bool one_sbi = this->mfds_parts.size() == 1 && this->mfds_parts[0].kind == "sbi" && observed.empty();
    if (reasons.empty() && likelihoods > 0 && !one_sbi && this->mfds_sbi_needs_direct)
      reasons.push_back("with several likelihood factors, each simulator must be a C++ function without summaries");
    std::vector<std::string> data_variables;
    std::string body;
    if (reasons.empty() && likelihoods > 0 && !one_sbi)
    {
      for (const std::string &y : observed)
      {
        std::set<std::string> read, nodes_seen;
        for (const ModelStatement &st : this->model_statements)
          if (st.op == "~" && st.variable == y)
          {
            this->data_read_by(st.value, read, nodes_seen);
            for (size_t i = 0; i < st.indices.size(); ++i) this->data_read_by(st.indices[i], read, nodes_seen);
            for (const ModelLoop &loop : st.loops) { this->data_read_by(loop.from, read, nodes_seen); this->data_read_by(loop.to, read, nodes_seen); }
          }
        for (const std::string &name : read)
          if (!this->is_literal_data(name))
            reasons.push_back("the statements of " + y + " use the data " + name +
                              ", which is given by a function (MFDS's simulator receives no data)");
      }
      for (const MfdsPart &part : this->mfds_parts)
      {
        Arguments arguments = this->split_arguments(part.call);
        if (part.kind == "sbi")
          continue;
        std::vector<std::string> order = (part.kind == "linear") ? std::vector<std::string>{"data", "state", "matrix", "covariance"}
                                                                 : std::vector<std::string>{"data", "function", "covariance"};
        for (size_t k = 0; k < arguments.positional.size() && k < order.size(); ++k)
          arguments.named[order[k]] = arguments.positional[k];
        if (mfds_names(arguments.named["data"]).size() != 1)
          reasons.push_back("the factor " + part.name + " has several data variables");
        if (this->mfds_matrix_code(arguments.named["covariance"]).empty() ||
            (part.kind == "linear" && this->mfds_matrix_code(arguments.named["matrix"]).empty()) ||
            (part.kind == "nonlinear" && (!this->is_function_reference(arguments.named["function"]) ||
                                         arguments.named["function"]["namespace"].as_string() != "cpp")))
          reasons.push_back("the factor " + part.name + " uses a function in another language");
      }
    }
    bool generated = false;
    if (reasons.empty() && likelihoods > 0 && !one_sbi)
    {
     // (a model MFDS cannot simulate must still compile: an error here is
     // MFDS's reason, not the file's)
     try
     {
      for (const std::string &y : observed)
      {
        std::string rows, cols;
        this->variable_extent(y, rows, cols);
        body += "output[" + quoted(y) + "] = arma::mat(" + rows + ", " + cols + ", arma::fill::zeros);\n";
        for (const ModelStatement &st : this->model_statements)
        {
          if (st.op != "~" || st.variable != y)
            continue;
          CodeContext context;
          context.prior = true;
          context.data = "output";
          for (const ModelLoop &loop : st.loops) context.indices.insert(loop.index);
          std::string assign = (st.value["function"].as_string() == "mvnorm")
                               ? "output[" + quoted(y) + "] = " + this->random_code(st, context) + ";\n"
                               : this->target_code(st, context, true) + " = " + this->random_code(st, context) + ";\n";
          CodeContext fresh;
          fresh.prior = true;
          body += this->loops_code(st, fresh, assign, "");
        }
        data_variables.push_back(y);
      }
      for (const MfdsPart &part : this->mfds_parts)
      {
        Arguments arguments = this->split_arguments(part.call);
        if (part.kind == "sbi")
        {
          std::string f = this->cpp_function(arguments.named["simulate"], "the simulator");
          body += "{\n  Data simulated = " + f + "(rng, parameters);\n";
          for (const std::string &v : mfds_names(arguments.named["data"]))
          {
            body += "  output[" + quoted(v) + "] = simulated[" + quoted(v) + "];\n";
            data_variables.push_back(v);
          }
          body += "}\n";
          continue;
        }
        std::vector<std::string> order = (part.kind == "linear") ? std::vector<std::string>{"data", "state", "matrix", "covariance"}
                                                                 : std::vector<std::string>{"data", "function", "covariance"};
        for (size_t k = 0; k < arguments.positional.size() && k < order.size(); ++k)
          arguments.named[order[k]] = arguments.positional[k];
        const std::string y = mfds_names(arguments.named["data"])[0]; // (const: the file's quoted, not std::quoted)
        body += "{\n";
        if (part.kind == "linear")
        {
          std::vector<std::string> state;
          for (const std::string &x : mfds_names(arguments.named["state"]))
            state.push_back("arma::vectorise(parameters[" + quoted(x) + "])");
          std::string stacked = state.empty() ? "arma::colvec()" : state[0];
          for (size_t k = 1; k < state.size(); ++k)
            stacked = "arma::join_cols(" + stacked + ", " + state[k] + ")";
          body += "  arma::colvec mean = arma::vectorise(" + this->mfds_matrix_code(arguments.named["matrix"]) + " * " + stacked + ");\n";
        }
        else
          body += "  arma::colvec mean = arma::vectorise(" + this->cpp_function(arguments.named["function"], "the function") +
                  "(parameters)[" + quoted(y) + "]);\n";
        body += "  arma::mat covariance = " + this->mfds_matrix_code(arguments.named["covariance"]) + ";\n"
                "  output[" + quoted(y) + "] = arma::mat(mean + rmvnorm(rng, arma::zeros<arma::colvec>(mean.n_elem), covariance));\n}\n";
        data_variables.push_back(y);
      }
      std::string indented, line;
      std::stringstream lines(body);
      while (std::getline(lines, line)) indented += "  " + line + "\n";
      this->mfds_code += "// MFDS: the likelihood's simulator (dev/plan.md 3.16)\n"
                         "Data ilike_mfds_simulate(RandomNumberGenerator &rng, const Parameters &parameters)\n{\n"
                         "  Data output;\n" + indented + "  return output;\n}\n\n";
      this->mfds_entry.push_back(this->placeholder("ilike_mfds_simulate", "mfds_simulate"), "simulate");
      this->mfds_entry.push_back(this->placeholder(this->generate_string_function(joined(data_variables, ";")), "data_variable"),
                                 "data_variable");
      generated = true;
     }
     catch (std::exception &e)
     {
       reasons.push_back(std::string("its simulator could not be generated (") + e.what() + ")");
     }
    }
    if (generated)
      ;
    else if (!reasons.empty())
      this->mfds_entry.push_back(RecipeValue(joined(reasons, "; ")), "simulate_reason");
    else if (likelihoods == 0)
      this->mfds_entry.push_back(RecipeValue(std::string("the model has no likelihood")), "simulate_reason");

    // the Gaussian prior's moments, from prior statements
    std::vector<std::string> log_variables = this->mfds_log_variables;
    std::vector<size_t> prior;
    for (size_t k = 0; k < this->model_statements.size(); ++k)
      if (this->model_statements[k].op == "~" && !this->model_statements[k].observed)
        prior.push_back(k);
    std::vector<std::string> prior_reasons = this->mfds_prior_reasons;
    std::vector<std::string> variables;
    std::map<std::string, std::string> family;
    for (size_t k : prior)
    {
      const ModelStatement &st = this->model_statements[k];
      std::string f = st.value["function"].as_string();
      if (f != "norm" && f != "mvnorm" && f != "lnorm")
      {
        prior_reasons.push_back("the prior of " + st.variable + " is " + f + ", not norm, mvnorm or lnorm");
        continue;
      }
      std::set<std::string> found;
      this->read_parameters(st.value, found);
      for (size_t i = 0; i < st.indices.size(); ++i) this->read_parameters(st.indices[i], found);
      for (const ModelLoop &loop : st.loops) { this->read_parameters(loop.from, found); this->read_parameters(loop.to, found); }
      for (const ModelLoop &loop : st.loops) found.erase(loop.index);
      found.erase(st.variable);
      if (!found.empty())
      {
        prior_reasons.push_back("the prior of " + st.variable + " depends on " + *found.begin() +
                                " (a prior with parameters in its arguments is not Gaussian)");
        continue;
      }
      std::string kind = (f == "lnorm") ? "lnorm" : "norm";
      if (family.count(st.variable) && family[st.variable] != kind)
        prior_reasons.push_back("the elements of " + st.variable + " mix log-normal and normal priors");
      if (!family.count(st.variable))
        variables.push_back(st.variable);
      family[st.variable] = kind;
    }
    bool prior_generated = false;
    if (!prior.empty() && prior_reasons.empty())
    {
     try
     {
      std::string code = "// MFDS: the Gaussian prior's means and covariances (dev/plan.md 3.16)\n"
                         "Data ilike_mfds_prior()\n{\n  Data output;\n"
                         "  Parameters parameters; // (for deterministic nodes: the prior's arguments read no parameters)\n";
      for (size_t v = 0; v < variables.size(); ++v)
      {
        const std::string &name = variables[v];
        bool whole = false;
        for (size_t k : prior)
        {
          const ModelStatement &st = this->model_statements[k];
          if (st.variable == name && st.value["function"].as_string() == "mvnorm")
          {
            CodeContext context;
            context.prior = true;
            const RecipeValue &arguments = st.value["arguments"];
            code += "  output[" + quoted(name) + "] = arma::mat(arma::vectorise(" + this->matrix_code(arguments[0], context) + "));\n"
                    "  output[" + quoted(name + ":cov") + "] = " + this->matrix_code(arguments[1], context) + ";\n";
            whole = true;
          }
        }
        if (!whole)
        {
          std::string rows, cols, mean = "mean_" + std::to_string(v), variance = "variance_" + std::to_string(v);
          this->variable_extent(name, rows, cols);
          code += "  arma::mat " + mean + "(" + rows + ", " + cols + ", arma::fill::zeros), " + variance + "(" + rows + ", " + cols +
                  ", arma::fill::zeros);\n";
          for (size_t k : prior)
          {
            const ModelStatement &st = this->model_statements[k];
            if (st.variable != name)
              continue;
            CodeContext context;
            context.prior = true;
            for (const ModelLoop &loop : st.loops) context.indices.insert(loop.index);
            std::string element = "(0)";
            if (st.indices.size() == 1)
              element = "(" + this->index_code(st.indices[0], context) + ")";
            else if (st.indices.size() == 2)
              element = "(" + this->index_code(st.indices[0], context) + ", " + this->index_code(st.indices[1], context) + ")";
            const RecipeValue &arguments = st.value["arguments"];
            std::string assign = mean + element + " = " + this->scalar_code(arguments[0], context) + ";\n" +
                                 variance + element + " = std::pow(" + this->scalar_code(arguments[1], context) + ", 2);\n";
            CodeContext fresh;
            fresh.prior = true;
            code += this->loops_code(st, fresh, assign, "  ");
          }
          code += "  if (arma::any(arma::vectorise(" + variance + ") <= 0.0))\n"
                  "    Rcpp::stop(" + quoted("MFDS: every element of " + name + " needs a prior with a positive standard deviation.") + ");\n"
                  "  output[" + quoted(name) + "] = " + mean + ";\n"
                  "  output[" + quoted(name + ":cov") + "] = arma::mat(arma::diagmat(arma::vectorise(" + variance + ")));\n";
        }
        if (family[name] == "lnorm")
        {
          code += "  output[" + quoted(name + ":log") + "] = arma::mat(1, 1, arma::fill::ones);\n";
          log_variables.push_back(name);
        }
      }
      code += "  return output;\n}\n\n";
      this->mfds_code += code;
      this->mfds_entry.push_back(this->placeholder("ilike_mfds_prior", "mfds_prior"), "prior");
      prior_generated = true;
     }
     catch (std::exception &e)
     {
       prior_reasons.push_back(std::string("its Gaussian moments could not be generated (") + e.what() + ")");
       log_variables = this->mfds_log_variables;
     }
    }
    if (prior_generated)
      ;
    else if (!prior_reasons.empty())
      this->mfds_entry.push_back(RecipeValue(joined(prior_reasons, "; ")), "prior_reason");

    // the map from the Gaussian scale to the model's, for log-normal variables
    if (!log_variables.empty())
    {
      std::string code = "// MFDS: from the Gaussian scale to the model's (dev/plan.md 3.16)\n"
                         "Parameters ilike_mfds_log_transform(const Parameters &u)\n{\n"
                         "  Parameters theta = u.deep_copy(); // (a plain copy shares u's matrices, which must not change)\n";
      for (const std::string &name : log_variables)
        code += "  theta[" + quoted(name) + "] = arma::exp(u[" + quoted(name) + "]);\n";
      code += "  return theta;\n}\n\n";
      this->mfds_code += code;
      this->mfds_entry.push_back(this->placeholder("ilike_mfds_log_transform", "mfds_transform"), "log_transform");
    }
  }

  // the size of a variable with "~" statements (a parameter, or data when
  // simulating it): the largest index of its statements, as C++ expressions
  void variable_extent(const std::string &variable, std::string &rows, std::string &cols)
  {
    rows = "1";
    cols = "1";
    for (const ModelStatement &other : this->model_statements)
    {
      if (other.op != "~" || other.variable != variable)
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
            this->error(other.node, "an index of " + variable + " in its statements is a loop index or a number (e.g. x[1], x[t]).");
          extent = std::to_string((int)number(0));
        }
        std::string &dimension = (i == 0) ? rows : cols;
        dimension = "std::max(" + dimension + ", " + extent + ")";
      }
    }
  }

  void finish_model_statements()
  {
    if (this->model_statements.empty())
      return;
    std::string code = R"ILIKE(// model statements (generated by the ilike language, I4a)
inline double ilike_model_dbeta(double x, double a, double b)
{
  if (x <= 0.0 || x >= 1.0) return -std::numeric_limits<double>::infinity();
  return std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b) + (a - 1.0)*std::log(x) + (b - 1.0)*std::log(1.0 - x);
}
// an element read by a model statement, checked: outside the variable, an
// index that reads parameters (a proposed value outside its range) makes the
// log density -Inf (ilike_model_outside, caught by the log density's
// function), and any other index is an error in the model
struct ilike_model_outside : public std::runtime_error
{
  ilike_model_outside(const std::string &what) : std::runtime_error(what) {}
};
// the model's file, for messages (set by R after loading: the generated code
// holds the marker @ILIKE_FILE@ instead, so that it does not depend on the path)
inline std::string &ilike_model_file()
{
  static std::string file;
  return file;
}
inline std::string ilike_model_where(const char *what)
{
  std::string text(what);
  size_t at = text.find("@ILIKE_FILE@");
  if (at != std::string::npos)
    text.replace(at, 12, ilike_model_file());
  return text;
}
inline void ilike_model_outside_index(const arma::mat &m, int i, int j, bool stochastic, const char *what)
{
  std::string message = ilike_model_where(what) + ": the index (" + std::to_string(i + 1) +
                        (j >= 0 ? ", " + std::to_string(j + 1) : std::string("")) + ") is outside the variable, which has " +
                        std::to_string(m.n_rows) + " row(s) and " + std::to_string(m.n_cols) + " column(s)";
  if (stochastic)
    throw ilike_model_outside(message + " (the index reads parameters, so this point has density 0).");
  Rcpp::stop(message + ".");
}
inline double ilike_model_at(const arma::mat &m, int i, bool stochastic, const char *what)
{
  if ((i < 0) || ((arma::uword)i >= m.n_elem))
    ilike_model_outside_index(m, i, -1, stochastic, what);
  return m(i);
}
inline double ilike_model_at(const arma::mat &m, int i, int j, bool stochastic, const char *what)
{
  if ((i < 0) || (j < 0) || ((arma::uword)i >= m.n_rows) || ((arma::uword)j >= m.n_cols))
    ilike_model_outside_index(m, i, j, stochastic, what);
  return m(i, j);
}
// gradients (I4e): an adjoint added to element i (or (i, j)) of the gradient
// of the variable asked for, with the same check as ilike_model_at
inline void ilike_model_add(arma::mat &g, int i, double a, bool stochastic, const char *what)
{
  if ((i < 0) || ((arma::uword)i >= g.n_elem))
    ilike_model_outside_index(g, i, -1, stochastic, what);
  g(i) += a;
}
inline void ilike_model_add(arma::mat &g, int i, int j, double a, bool stochastic, const char *what)
{
  if ((i < 0) || (j < 0) || ((arma::uword)i >= g.n_rows) || ((arma::uword)j >= g.n_cols))
    ilike_model_outside_index(g, i, j, stochastic, what);
  g(i, j) += a;
}
// the variables whose gradients are asked for, and their gradients (C204):
// find(name) is the gradient of name, or NULL if it is not asked for
struct ilike_model_gradients
{
  const std::vector<std::string> &names;
  std::vector<arma::mat> &values;
  std::unordered_map<std::string, size_t> index; // (only for many variables)
  ilike_model_gradients(const std::vector<std::string> &names_in, std::vector<arma::mat> &values_in,
                        const Parameters &parameters)
  : names(names_in), values(values_in)
  {
    this->values.resize(this->names.size());
    for (size_t k = 0; k < this->names.size(); ++k)
      this->values[k].zeros(arma::size(parameters[this->names[k]]));
    if (this->names.size() > 8)
      for (size_t k = 0; k < this->names.size(); ++k)
        this->index[this->names[k]] = k;
  }
  arma::mat* find(const char *name)
  {
    if (this->index.empty())
    {
      for (size_t k = 0; k < this->names.size(); ++k)
        if (this->names[k] == name)
          return &this->values[k];
      return NULL;
    }
    auto found = this->index.find(name);
    return (found == this->index.end()) ? NULL : &this->values[found->second];
  }
};
// a phi(a) and phi(a) for the truncated normal, 0 at infinite limits
inline double ilike_model_phi(double a)
{
  return std::isfinite(a) ? R::dnorm(a, 0.0, 1.0, 0) : 0.0;
}
inline double ilike_model_a_phi(double a)
{
  return std::isfinite(a) ? a*R::dnorm(a, 0.0, 1.0, 0) : 0.0;
}
// the partial derivatives of log dtnorm(x, m, s, l, u) with respect to x, m, s, l, u
inline arma::colvec ilike_model_dtnorm_partials(double x, double m, double s, double l, double u)
{
  double alpha = (l - m)/s, beta = (u - m)/s;
  double Z = R::pnorm(beta, 0.0, 1.0, 1, 0) - R::pnorm(alpha, 0.0, 1.0, 1, 0);
  arma::colvec d(5);
  d[0] = -(x - m)/(s*s);
  d[1] = (x - m)/(s*s) + (ilike_model_phi(beta) - ilike_model_phi(alpha))/(s*Z);
  d[2] = -1.0/s + (x - m)*(x - m)/(s*s*s) + (ilike_model_a_phi(beta) - ilike_model_a_phi(alpha))/(s*Z);
  d[3] = ilike_model_phi(alpha)/(s*Z);
  d[4] = -ilike_model_phi(beta)/(s*Z);
  return d;
}
inline double ilike_model_rbeta(RandomNumberGenerator &rng, double a, double b)
{
  double x = rgamma(rng, a, 1.0), y = rgamma(rng, b, 1.0);
  return x/(x + y);
}
inline double ilike_model_rbern(RandomNumberGenerator &rng, double p)
{
  return (runif(rng) < p) ? 1.0 : 0.0;
}
inline double ilike_model_rbinom(RandomNumberGenerator &rng, double n, double p)
{
  double count = 0.0;
  for (int k = 0; k < (int)n; ++k)
    if (runif(rng) < p) count = count + 1.0;
  return count;
}
inline double ilike_model_dbern(double x, double p)
{
  if (x == 1.0) return std::log(p);
  if (x == 0.0) return std::log(1.0 - p);
  return -std::numeric_limits<double>::infinity();
}
inline void ilike_model_read(std::vector<double> &reads, double term, double variable, double row, double col)
{
  reads.push_back(term); reads.push_back(variable); reads.push_back(row); reads.push_back(col);
}
inline arma::mat ilike_model_table(const std::vector<double> &values, size_t width)
{
  arma::mat table(values.size()/width, width);
  for (size_t r = 0; r < table.n_rows; ++r)
    for (size_t c = 0; c < width; ++c)
      table(r, c) = values[r*width + c];
  return table;
}
inline double ilike_model_dbinom(double x, double n, double p)
{
  if (x < 0.0 || x > n || x != std::floor(x)) return -std::numeric_limits<double>::infinity();
  return std::lgamma(n + 1.0) - std::lgamma(x + 1.0) - std::lgamma(n - x + 1.0) + x*std::log(p) + (n - x)*std::log(1.0 - p);
}

)ILIKE";
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
      context.node = true;
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
      // what the node reads (I4b), with the same arguments
      std::string reads_signature = "inline void ilike_model_node_reads_" + st.variable +
                                    "(std::vector<double> &reads, double term, const Data *data";
      for (size_t k = 0; k < st.indices.size(); ++k)
        reads_signature += ", int " + st.indices[k]["name"].as_string();
      node_code += reads_signature + ");\n";
      this->deferred_node_code += reads_signature + ")\n{\n  Parameters parameters; // (no index read here depends on the parameters)\n" +
                                  this->reads_code(st.value, context, false) + "}\n\n";
      // its gradient (I4e): the adjoint passed on to what the node reads
      std::string gradient_signature = "inline void ilike_model_node_gradient_" + st.variable +
                                       "(ilike_model_gradients &ilike_G, double ilike_adjoint, "
                                       "const Parameters &parameters, const Data *data";
      for (size_t k = 0; k < st.indices.size(); ++k)
        gradient_signature += ", int " + st.indices[k]["name"].as_string();
      gradient_signature += ")";
      node_code += gradient_signature + ";\n";
      std::string gradient_body;
      try
      {
        gradient_body = this->gradient_code(st.value, "ilike_adjoint", context);
      }
      catch (std::exception &e)
      {
        gradient_body = "Rcpp::stop(" + quoted("The gradient through the deterministic node " + st.variable + " is not available.") + ");\n";
      }
      this->deferred_node_code += gradient_signature + "\n{\n" + indented_code(gradient_body) + "}\n\n";
    }
    code += node_code + "\n";

    // a factor's name is not a parameter
    for (const ModelStatement &st : this->model_statements)
      if (st.op == "+=")
        for (const ModelStatement &other : this->model_statements)
          if (other.op != "+=" && other.variable == st.variable)
            this->error(st.node, st.variable + " is a variable in the model statements, so it cannot name a factor of terms.");

    // the likelihood factors: one per observed variable, and one per name of
    // "factor name += expression" terms
    struct PendingPlate { size_t position; std::string suffix; std::string name; std::string gradient; };
    std::vector<PendingPlate> plates;
    std::vector<std::string> observed_done;
    for (const ModelStatement &first : this->model_statements)
    {
      if (!this->likelihood_term(first) ||
          std::find(observed_done.begin(), observed_done.end(), first.variable) != observed_done.end())
        continue;
      observed_done.push_back(first.variable);
      std::string function = "ilike_model_log_likelihood_" + first.variable;
      std::string body = "double log_density = 0.0;\n";
      for (const ModelStatement &st : this->model_statements)
      {
        if (!this->likelihood_term(st) || st.variable != first.variable)
          continue;
        // (the loop indices in scope for the target and the arguments)
        CodeContext inner_context;
        for (const ModelLoop &loop : st.loops) inner_context.indices.insert(loop.index);
        std::string statement_body = "log_density += " + this->term_code(st, inner_context) + ";\n";
        CodeContext fresh;
        body += this->loops_code(st, fresh, statement_body, "");
      }
      body += "return log_density;\n";
      std::string indented;
      std::stringstream lines(body);
      std::string line;
      while (std::getline(lines, line)) indented += "  " + line + "\n";
      this->deferred_function_code += "double " + function + "(const Parameters &parameters, const Data &data)\n{\n  try\n  {\n" +
                                      indented + "  }\n" + std::string("  catch (const ilike_model_outside &)\n  {\n    return -std::numeric_limits<double>::infinity(); // (an index that reads parameters is outside its variable)\n  }\n") + "}\n\n";
      RecipeValue factor = RecipeValue::list();
      factor.push_back(this->placeholder(function, "evaluate_log_likelihood"), "evaluate_log_likelihood");

      // its plate (I4b): the same terms, in the same order
      std::vector<const ModelStatement*> statements;
      for (const ModelStatement &st : this->model_statements)
        if (this->likelihood_term(st) && st.variable == first.variable)
          statements.push_back(&st);
      // and its gradient (I4e)
      std::string gradient = this->generate_gradient("likelihood_" + first.variable, statements, false);
      if (!gradient.empty())
        factor.push_back(this->placeholder(gradient, "evaluate_gradient_log_likelihood"), "evaluate_gradient_log_likelihood");
      this->factor_entries[this->model_factor_position.at(first.variable)] = factor;
      this->generated_factor_reads[this->factor_numbers.at(first.variable)] = this->statements_read(statements);
      if (this->plate_possible(statements))
      {
        this->generate_plate("likelihood_" + first.variable, statements, false);
        plates.push_back({this->model_factor_position.at(first.variable), "likelihood_" + first.variable, first.variable, gradient});
      }
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
            std::string rows, cols;
            this->variable_extent(st.variable, rows, cols);
            simulate_body += "parameters[" + quoted(st.variable) + "] = arma::mat(" + rows + ", " + cols + ", arma::fill::zeros);\n";
          }
        }
        std::string assign;
        if (st.value["function"].as_string() == "mvnorm")
          assign = "parameters[" + quoted(st.variable) + "] = " + this->random_code(st, loop_context) + ";\n";
        else
          assign = this->target_code(st, loop_context, true) + " = " + this->random_code(st, loop_context) + ";\n";
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
      this->deferred_function_code += "double ilike_model_log_prior(const Parameters &parameters)\n{\n  try\n  {\n" + indent(log_body) +
                                      "  }\n" + std::string("  catch (const ilike_model_outside &)\n  {\n    return -std::numeric_limits<double>::infinity(); // (an index that reads parameters is outside its variable)\n  }\n") + "}\n\n";
      this->deferred_function_code += "Parameters ilike_model_simulate_prior(RandomNumberGenerator &rng)\n{\n" + indent(simulate_body) + "}\n\n";
      RecipeValue factor = RecipeValue::list();
      factor.push_back(this->placeholder("ilike_model_log_prior", "evaluate_log_prior"), "evaluate_log_prior");
      factor.push_back(this->placeholder("ilike_model_simulate_prior", "simulate_prior"), "simulate_prior");

      std::vector<const ModelStatement*> statements;
      for (size_t k : order)
        statements.push_back(&this->model_statements[k]);
      // its gradient (I4e)
      std::string gradient = this->generate_gradient("prior", statements, true);
      if (!gradient.empty())
        factor.push_back(this->placeholder(gradient, "evaluate_gradient_log_prior"), "evaluate_gradient_log_prior");
      this->factor_entries[this->model_factor_position.at("model_prior")] = factor;
      this->generated_factor_reads[this->factor_numbers.at("model_prior")] = this->statements_read(statements);
      if (this->plate_possible(statements))
      {
        this->generate_plate("prior", statements, true);
        plates.push_back({this->model_factor_position.at("model_prior"), "prior", "model_prior",
                          gradient.empty() ? gradient : gradient + "_plate"});
      }
    }

    // the plates' entries, now that every parameter read has its number
    for (const PendingPlate &plate : plates)
    {
      RecipeValue entry = this->plate_entry(plate.suffix, plate.name);
      if (!plate.gradient.empty())
      {
        entry.push_back(this->placeholder(plate.gradient, "evaluate_gradient_log_likelihood"), "gradient");
        // the gradients of several variables in one pass (C204): f_many, or
        // for the prior f_many_plate (f_plate's counterpart)
        std::string many = plate.gradient;
        if (many.size() > 6 && many.compare(many.size() - 6, 6, "_plate") == 0)
          many = many.substr(0, many.size() - 6) + "_many_plate";
        else
          many = many + "_many";
        entry.push_back(this->placeholder(many, "plate_gradients"), "gradients");
      }
      RecipeValue variables = RecipeValue::list();
      for (const std::string &name : this->plate_variable_names)
        variables.push_back(RecipeValue(name));
      entry.push_back(variables, "variables");
      this->factor_entries[plate.position].push_back(entry, "plate");
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
    // (without a method the factor is a simulator alone, used by MFDS)
    if (!method.is_null() && !role_exists(this->builtin_name(method), BuiltinRole::SBILikelihood))
      find_builtin(this->builtin_name(method), BuiltinRole::SBILikelihood); // (the registry's message)
    if (!arguments.has("simulate"))
      this->error(call, "sbi: the simulator is missing (simulate = cpp::f).");
    if (!arguments.has("data"))
      this->error(call, "sbi: the data variable is missing (data = y).");
    RecipeValue factor = RecipeValue::list();
    factor.push_back(this->function_placeholder(arguments.named["simulate"], "simulate_data_model", "simulate"),
                     "simulate_data_model");
    if (!method.is_null())
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

  // a copy of an expression with the loop indices replaced by their values
  static RecipeValue substitute(const RecipeValue &expression, const std::map<std::string, double> &values)
  {
    if (!expression.is_list())
      return expression;
    if (kind_of(expression) == "name" && expression["namespace"].as_string().empty())
    {
      auto found = values.find(expression["name"].as_string());
      if (found != values.end())
      {
        RecipeValue number = RecipeValue::list();
        number.push_back(RecipeValue("number"), "expression");
        number.push_back(RecipeValue(found->second), "value");
        if (expression.contains("line")) number.push_back(expression["line"], "line");
        if (expression.contains("column")) number.push_back(expression["column"], "column");
        return number;
      }
    }
    RecipeValue copy = RecipeValue::list();
    std::vector<std::string> names = expression.names();
    for (size_t k = 0; k < expression.size(); ++k)
      copy.push_back(substitute(expression[k], values), (k < names.size()) ? names[k] : std::string());
    return copy;
  }

  // for (j in a:b) move name = metropolis(norm_rw(theta[j], s)) (I4c): a scan
  // over the elements, in the loops' order (the last loop fastest), as one
  // move whose entry lists the elements (rows, columns, scales; 0-based)
  // the constants named anywhere in an expression or statement
  void constants_named(const RecipeValue &node, std::set<std::string> &found) const
  {
    if (!node.is_list())
      return;
    if (kind_of(node) == "name" && node["namespace"].as_string().empty() &&
        this->constant_names.count(node["name"].as_string()))
      found.insert(node["name"].as_string());
    for (size_t k = 0; k < node.size(); ++k)
      this->constants_named(node[k], found);
  }

  // --- Gibbs updates from normal full conditionals (C206) ---------------------

  bool reads_name(const RecipeValue &expression, const std::string &x) const
  {
    std::set<std::string> found;
    this->read_parameters(expression, found);
    return found.count(x) > 0;
  }

  // whether an expression is affine in x (a + b x, with a and b not depending
  // on x), looking through deterministic nodes; `why` says why not
  bool affine_in(const RecipeValue &expression, const std::string &x, std::string &why) const
  {
    if (!this->reads_name(expression, x))
      return true;
    std::string kind = kind_of(expression);
    if (kind == "name" || kind == "index")
    {
      std::string name = (kind == "name") ? expression["name"].as_string() : expression["target"]["name"].as_string();
      if (kind == "index")
        for (size_t k = 0; k < expression["indices"].size(); ++k)
          if (this->reads_name(expression["indices"][k], x))
          {
            why = "an index depends on " + x;
            return false;
          }
      if ((name != x) && this->is_deterministic(name))
        return this->affine_in(this->deterministic_statement(name).value, x, why);
      return true;
    }
    if (kind == "unary")
      return this->affine_in(expression["operand"], x, why);
    if (kind == "binary")
    {
      std::string op = expression["operator"].as_string();
      const RecipeValue &left = expression["left"], &right = expression["right"];
      if ((op == "+") || (op == "-"))
        return this->affine_in(left, x, why) && this->affine_in(right, x, why);
      if (op == "*")
      {
        if (!this->reads_name(left, x))
          return this->affine_in(right, x, why);
        if (!this->reads_name(right, x))
          return this->affine_in(left, x, why);
        why = "two factors of a product depend on " + x;
        return false;
      }
      if (op == "/")
      {
        if (this->reads_name(right, x))
        {
          why = "it divides by an expression depending on " + x;
          return false;
        }
        return this->affine_in(left, x, why);
      }
      why = x + " is in a power (^)";
      return false;
    }
    if (kind == "call")
    {
      why = x + " is inside " + expression["function"].as_string() + "(…)";
      return false;
    }
    why = "it is not a sum of multiples of " + x;
    return false;
  }

  // Whether x's full conditional is normal by the model statements: x's own
  // distribution is norm or mvnorm with arguments not depending on x, and
  // every other term depending on x is a norm whose mean is affine in x and
  // whose sd, target and loops do not depend on x. Then its log density is
  // exactly quadratic in each element of x. The answer: "" when it is, or why
  // not (with the statement at fault, if any).
  std::pair<const ModelStatement*, std::string> normal_conditional_problem(const std::string &x) const
  {
    bool own = false;
    for (const ModelStatement &st : this->model_statements)
    {
      if (st.op == "<-")
        continue;
      if (!this->statements_read(std::vector<const ModelStatement*>(1, &st)).count(x))
        continue;
      if (st.op == "+=")
        return {&st, "the term of factor " + st.variable + " depends on " + x};
      std::string f = st.value["function"].as_string();
      const RecipeValue &arguments = st.value["arguments"];
      for (const ModelLoop &loop : st.loops)
        if (this->reads_name(loop.from, x) || this->reads_name(loop.to, x))
          return {&st, "a loop's range depends on " + x};
      for (size_t k = 0; k < st.indices.size(); ++k)
        if (this->reads_name(st.indices[k], x))
          return {&st, "an index depends on " + x};
      if (!st.observed && (st.variable == x))
      {
        own = true;
        if ((f != "norm") && (f != "mvnorm"))
          return {&st, x + "'s distribution is " + f + " (a normal full conditional needs norm or mvnorm)"};
        for (size_t k = 0; k < arguments.size(); ++k)
          if (this->reads_name(arguments[k], x))
            return {&st, "the arguments of " + x + "'s distribution depend on " + x};
        continue;
      }
      if (f != "norm")
        return {&st, "this term's distribution is " + f + " (a normal full conditional needs norm terms)"};
      if (this->reads_name(arguments[1], x))
        return {&st, x + " is in the standard deviation"};
      std::string why;
      if (!this->affine_in(arguments[0], x, why))
        return {&st, x + " enters the mean nonlinearly: " + why};
    }
    if (!own)
      return {NULL, x + " has no distribution in the model statements"};
    return {NULL, ""};
  }

  Move element_move_loop(const RecipeValue &statement)
  {
    // (its constants take their values now, when compiling)
    this->constants_named(statement, this->element_constants_used);
    // the loops, down to the one move they hold
    std::vector<ModelLoop> loops;
    const RecipeValue *current = &statement;
    while (current->contains("keyword") && (*current)["keyword"].as_string() == "for")
    {
      const RecipeValue &range = (*current)["range"];
      if (kind_of(range) != "binary" || range["operator"].as_string() != ":")
        this->error(*current, "a loop's range is from:to (e.g. for (j in 1:J)).");
      ModelLoop loop;
      loop.index = (*current)["index"].as_string();
      loop.from = range["left"];
      loop.to = range["right"];
      loops.push_back(loop);
      const RecipeValue &body = (*current)["body"];
      if (body.size() != 1)
        this->error(*current, "a loop in an algorithm holds one move (for (j in 1:J) move m = metropolis(norm_rw(theta[j], s))).");
      current = &body[0];
    }
    const RecipeValue &move_statement = *current;
    if (!(move_statement.contains("keyword") && move_statement["keyword"].as_string() == "move") ||
        move_statement["operator"].as_string() != "=")
      this->error(move_statement, "a loop in an algorithm holds one move (move name = metropolis(…)).");
    std::string name = move_statement["target"]["name"].as_string();
    const RecipeValue &call = move_statement["value"];
    if (kind_of(call) != "call" || (call["function"].as_string() != "metropolis" && call["function"].as_string() != "gibbs"))
      this->error(move_statement, "a move on elements is a metropolis move (norm_rw(theta[j], sd), unif_rw(theta[j], halfwidth) "
                  "or flip(s[j])) or a gibbs move (gibbs(s[j])).");
    bool gibbs = (call["function"].as_string() == "gibbs");
    Arguments arguments = this->split_arguments(call);
    if (arguments.positional.empty())
      this->error(call, call["function"].as_string() + ": the element is missing.");
    std::string type, scale_name;
    RecipeValue element, scale, support;
    if (gibbs)
    {
      if (arguments.has("adapt") || arguments.has("smc_adapt"))
        this->error(call, "gibbs moves have no step to adapt.");
      this->allow_only(call, arguments, {"values", "factors", "cut_factors"}, 1);
      type = "gibbs";
      element = arguments.positional[0];
      if (arguments.has("values"))
        support = arguments.named["values"];
    }
    else
    {
      const RecipeValue &proposal = arguments.positional[0];
      type = (kind_of(proposal) == "call") ? this->builtin_name(proposal) : std::string();
      if (type != "norm_rw" && type != "unif_rw" && type != "flip")
        this->error(proposal, "a metropolis move on elements uses norm_rw(theta[j], sd), unif_rw(theta[j], halfwidth) or flip(s[j]).");
      Arguments proposal_arguments = this->split_arguments(proposal);
      scale_name = (type == "norm_rw") ? "sd" : "halfwidth";
      element = !proposal_arguments.positional.empty() ? proposal_arguments.positional[0]
              : (proposal_arguments.has("variable") ? proposal_arguments.named["variable"] : RecipeValue());
      if (type == "flip")
      {
        if (proposal_arguments.positional.size() > 1 || !proposal_arguments.named.empty())
          this->error(proposal, "flip(s[j]) takes only the element.");
      }
      else
      {
        scale = (proposal_arguments.positional.size() > 1) ? proposal_arguments.positional[1]
              : (proposal_arguments.has(scale_name) ? proposal_arguments.named[scale_name] : RecipeValue());
        if (scale.is_null())
          this->error(proposal, type + ": the " + scale_name + " is missing.");
      }
    }
    if (gibbs && (kind_of(element) != "name") && (kind_of(element) != "index"))
      this->error(call, "gibbs takes a scalar parameter (move m = gibbs(mu)) or an element in a loop "
                  "(for (j in 1:J) move m = gibbs(theta[j])); a proposal such as norm_rw(…) goes in metropolis(…).");
    bool whole_scalar = gibbs && loops.empty() && (kind_of(element) == "name");
    if (whole_scalar)
    {
      // gibbs(mu) on a whole variable: a scalar, as its element mu[1] (C206)
      std::string name = element["name"].as_string();
      for (const ModelStatement &st : this->model_statements)
        if ((st.op == "~") && !st.observed && (st.variable == name) &&
            ((st.indices.size() > 0) || st.value["function"].as_string() == "mvnorm"))
          this->error(call, "gibbs(" + name + ") moves a scalar; for the elements of " + name +
                      ", write for (j in 1:J) move m = gibbs(" + name + "[j]).");
      RecipeValue one = RecipeValue::list();
      one.push_back(RecipeValue("number"), "expression");
      one.push_back(RecipeValue(1.0), "value");
      RecipeValue indices = RecipeValue::list();
      indices.push_back(one, "");
      RecipeValue indexed = RecipeValue::list();
      indexed.push_back(RecipeValue("index"), "expression");
      indexed.push_back(element, "target");
      indexed.push_back(indices, "indices");
      element = indexed;
    }
    if (element.is_null() || kind_of(element) != "index" || kind_of(element["target"]) != "name")
      this->error(call, type + ": the element is written theta[j] (or x[i, j]) with the loop's indices.");
    std::string variable = element["target"]["name"].as_string();
    const RecipeValue &indices = element["indices"];
    if (indices.size() > 2)
      this->error(element, "at most two indices.");
    if (this->data_names.count(variable) || this->constant_names.count(variable))
      this->error(element, variable + " is not a parameter.");

    // every combination of the loops' values (each range evaluated with the
    // outer loops' values, so a range may depend on an outer index)
    std::vector<double> rows, columns, scales;
    // (scalar constants by their values when compiling)
    std::map<std::string, double> values;
    for (const std::string &constant : this->constant_names)
    {
      double number;
      if (this->constant_number(constant, number))
        values[constant] = number;
    }
    std::function<void(size_t)> enumerate = [&](size_t depth) {
      if (depth == loops.size())
      {
        arma::mat number;
        double index[2] = {1.0, 1.0};
        for (size_t k = 0; k < indices.size(); ++k)
        {
          if (!this->numeric(substitute(indices[k], values), number) || number.n_elem != 1)
            this->error(indices[k], "an element's index must be a number given the loops' indices.");
          index[k] = number(0);
        }
        if (index[0] < 1 || index[1] < 1 || index[0] != std::floor(index[0]) || index[1] != std::floor(index[1]))
          this->error(element, "an element's indices are whole numbers from 1.");
        double step = 0.0;
        if (!scale.is_null())
        {
          if (!this->numeric(substitute(scale, values), number) || number.n_elem != 1)
            this->error(scale, "the " + scale_name + " must be a number (it may use the loops' indices and constants).");
          step = number(0);
        }
        rows.push_back(index[0] - 1);
        columns.push_back(index[1] - 1);
        scales.push_back(step);
        return;
      }
      const ModelLoop &loop = loops[depth];
      arma::mat from, to;
      if (!this->numeric(substitute(loop.from, values), from) || !this->numeric(substitute(loop.to, values), to) ||
          from.n_elem != 1 || to.n_elem != 1)
        this->error(statement, "a loop's range in an algorithm must be numbers (constants allowed).");
      for (double v = from(0); v <= to(0); v += 1.0)
      {
        values[loop.index] = v;
        enumerate(depth + 1);
      }
      values.erase(loop.index);
    };
    enumerate(0);
    if (rows.empty())
      this->error(statement, "the loop has no elements.");

    // gibbs: the element's states (values = [...], or from the variable's
    // prior: bern, or binom with its n known when compiling)
    std::vector<double> states;
    if (gibbs)
    {
      arma::mat number;
      std::map<std::string, double> constants;
      for (const std::string &constant : this->constant_names)
      {
        double value;
        if (this->constant_number(constant, value))
          constants[constant] = value;
      }
      if (!support.is_null())
      {
        if (!this->numeric(substitute(support, constants), number) || number.n_elem == 0)
          this->error(support, "gibbs: values = [...] lists the element's possible values (numbers).");
        states.assign(number.begin(), number.end());
      }
      else
      {
        for (const ModelStatement &st : this->model_statements)
          if (st.op == "~" && !st.observed && st.variable == variable)
          {
            std::string f = st.value["function"].as_string();
            if (f == "bern")
              states = {0.0, 1.0};
            else if (f == "binom" && this->numeric(substitute(st.value["arguments"][0], constants), number) &&
                     number.n_elem == 1)
              for (double v = 0.0; v <= number(0); v += 1.0)
                states.push_back(v);
            break;
          }
        if (states.empty())
        {
          // not discrete: a normal full conditional, if the model statements
          // make it one (C206)
          std::pair<const ModelStatement*, std::string> problem = this->normal_conditional_problem(variable);
          if (!problem.second.empty())
          {
            std::string message = "gibbs(" + variable + (whole_scalar ? "" : "[…]") + "): " + variable + " is not discrete (bern, binom, or values = [...]), "
                                  "and its full conditional is not normal: " + problem.second +
                                  ". Use a metropolis move instead.";
            if (problem.first != NULL)
              this->error(problem.first->node, message);
            this->error(call, message);
          }
          // every factor reading it is known (what other factors read is not)
          if (!arguments.has("factors"))
          {
            std::vector<size_t> numbers = this->automatic_factor_index(variable);
            size_t total = this->factor_entries.size();
            for (size_t n = 1; n <= total; ++n)
              if ((numbers.empty() || std::find(numbers.begin(), numbers.end(), n) != numbers.end()) &&
                  !this->generated_factor_reads.count(n))
                this->error(call, "gibbs(" + variable + (whole_scalar ? "" : "[…]") + "): factor " + std::to_string(n) + " is not written as model "
                            "statements, so whether its density is normal in " + variable + " is not known; if it does not "
                            "depend on " + variable + ", list the factors the move uses with factors = [...].");
          }
          type = "gibbs_normal";
        }
      }
    }

    // the move's entry: as a metropolis move on the whole variable (which
    // gives its factors, from what the factors read, or from factors = [...]),
    // with the elements listed for the engine
    auto number_node = [](double value) {
      RecipeValue node = RecipeValue::list();
      node.push_back(RecipeValue("number"), "expression");
      node.push_back(RecipeValue(value), "value");
      return node;
    };
    auto call_node = [&](const std::string &function, const RecipeValue &call_arguments) {
      RecipeValue node = RecipeValue::list();
      node.push_back(RecipeValue("call"), "expression");
      node.push_back(RecipeValue(function), "function");
      node.push_back(RecipeValue(""), "namespace");
      node.push_back(call_arguments, "arguments");
      if (call.contains("line")) node.push_back(call["line"], "line");
      if (call.contains("column")) node.push_back(call["column"], "column");
      return node;
    };
    RecipeValue walk_arguments = RecipeValue::list();
    walk_arguments.push_back(element["target"], "");
    walk_arguments.push_back(number_node(scale.is_null() ? 1.0 : scales[0]), "");
    RecipeValue move_arguments = RecipeValue::list();
    move_arguments.push_back(call_node("norm_rw", walk_arguments), "");
    for (const auto &named : arguments.named)
      if (named.first == "factors" || named.first == "cut_factors")
        move_arguments.push_back(named.second, named.first);
    RecipeValue whole_call = call_node("metropolis", move_arguments);
    Move move = this->lower_move(move_statement, name, whole_call);
    RecipeValue elements = RecipeValue::list();
    elements.push_back(RecipeValue(variable), "variable");
    elements.push_back(RecipeValue(type), "type");
    elements.push_back(RecipeValue(rows), "rows");
    elements.push_back(RecipeValue(columns), "columns");
    elements.push_back(RecipeValue(scales), "scales");
    if (gibbs && (type == "gibbs"))
      elements.push_back(RecipeValue(states), "values");
    // adaptation of each element's step (dev/plan.md, "adaptation of element moves")
    if (arguments.has("adapt") || arguments.has("smc_adapt"))
    {
      if (type == "flip")
        this->error(call, "flip moves have no step to adapt.");
      if (arguments.has("adapt"))
      {
        const RecipeValue &adaptor = arguments.named["adapt"];
        if (kind_of(adaptor) != "call" || this->builtin_name(adaptor) != "scale")
          this->error(adaptor, "moves on elements adapt each element's step with adapt = scale(target, kappa, C, averaged).");
        RecipeValue until;
        RecipeValue scale_adaptor = this->without_until(adaptor, until);
        elements.push_back(this->builtin_node(scale_adaptor, true, variable), "adapt");
        if (!until.is_null())
          elements.push_back(until, "adapt_until");
      }
      if (arguments.has("smc_adapt"))
      {
        const RecipeValue &adaptor = arguments.named["smc_adapt"];
        if (kind_of(adaptor) != "call" || this->builtin_name(adaptor) != "gaussian")
          this->error(adaptor, "moves on elements are set from SMC particles with smc_adapt = gaussian(c).");
        elements.push_back(this->builtin_node(adaptor, true, variable), "smc_adapt");
      }
    }
    move.entry.push_back(elements, "m_elements");
    return move;
  }

  // --- transforms of moved variables (dev/plan.md 3.18 H2, C209) ---------------

  // a number known when compiling (constants by their values), or FALSE
  bool compile_time_number(const RecipeValue &expression, double &value)
  {
    std::map<std::string, double> constants;
    for (const std::string &constant : this->constant_names)
    {
      double c;
      if (this->constant_number(constant, c))
        constants[constant] = c;
    }
    arma::mat number;
    if (!this->numeric(substitute(expression, constants), number) || (number.n_elem != 1))
      return false;
    value = number(0);
    return true;
  }

  static RecipeValue transform_item(const std::string &variable, const std::string &kind, double lower, double upper)
  {
    RecipeValue item = RecipeValue::list();
    item.push_back(RecipeValue(variable), "variable");
    item.push_back(RecipeValue(kind), "kind");
    item.push_back(RecipeValue(lower), "lower");
    item.push_back(RecipeValue(upper), "upper");
    return item;
  }

  // the distribution of a parameter in the model section: its call
  // (f(arguments)), from the model statements, or from a simple prior
  // (theta ~ gamma(2, 1) with constant arguments, which stays the registry's
  // prior when no statement needs generated code); NULL if none
  RecipeValue distribution_of(const std::string &variable) const
  {
    for (const ModelStatement &st : this->model_statements)
      if ((st.op == "~") && !st.observed && (st.variable == variable))
        return st.value;
    for (const std::string &name : this->model_parameter_names())
      if (name == variable)
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
              if (st["keyword"].as_string().empty() && st.contains("operator") && (st["operator"].as_string() == "~") &&
                  (kind_of(st["target"]) == "name") && (st["target"]["name"].as_string() == variable) &&
                  (kind_of(st["value"]) == "call"))
                return st["value"];
            }
          }
        }
    return RecipeValue();
  }

  // the parameters given a distribution in the model section (statements,
  // and simple priors written "theta ~ f(…)")
  std::vector<std::string> model_parameter_names() const
  {
    std::vector<std::string> names;
    auto add = [&](const std::string &name) {
      if (!this->data_names.count(name) && (std::find(names.begin(), names.end(), name) == names.end()))
        names.push_back(name);
    };
    for (const ModelStatement &st : this->model_statements)
      if ((st.op == "~") && !st.observed)
        add(st.variable);
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
          if (st["keyword"].as_string().empty() && st.contains("operator") && (st["operator"].as_string() == "~") &&
              (kind_of(st["target"]) == "name"))
            add(st["target"]["name"].as_string());
        }
      }
    }
    return names;
  }

  // the transform of a variable from its support in the model section
  // (transform = auto): NULL when it needs none (or the support is not known)
  RecipeValue automatic_transform(const std::string &variable)
  {
    RecipeValue call = this->distribution_of(variable);
    if (!call.is_null())
    {
      std::string f = call["function"].as_string();
      const RecipeValue &arguments = call["arguments"];
      if ((f == "gamma") || (f == "lnorm") || (f == "exp"))
        return transform_item(variable, "log", 0.0, 0.0);
      if (f == "beta")
        return transform_item(variable, "logit", 0.0, 1.0);
      double a, b;
      if ((f == "unif") && this->compile_time_number(arguments[0], a) && this->compile_time_number(arguments[1], b))
        return transform_item(variable, "logit", a, b);
      if ((f == "tnorm") && (arguments.size() == 4) && this->compile_time_number(arguments[2], a) &&
          this->compile_time_number(arguments[3], b))
      {
        if (std::isfinite(a) && std::isfinite(b))
          return transform_item(variable, "logit", a, b);
        if (std::isfinite(a))
          return transform_item(variable, "lower", a, 0.0);
        if (std::isfinite(b))
          return transform_item(variable, "upper", 0.0, b);
      }
    }
    return RecipeValue();
  }

  // transform = auto | none | [log(x), logit(p), logit(p, a, b), lower(x, a), upper(x, b)]
  // for the moved variables (a stacked name's parts); NULL for none
  RecipeValue move_transform(const RecipeValue *argument, const std::string &variables, bool automatic_by_default)
  {
    // (the stacked name's parts, "a;b;c")
    std::vector<std::string> moved;
    std::stringstream parts(variables);
    for (std::string part; std::getline(parts, part, ';');)
      moved.push_back(part);
    bool automatic = automatic_by_default;
    RecipeValue list = RecipeValue::list();
    if (argument != NULL)
    {
      const RecipeValue &given = *argument;
      if ((kind_of(given) == "name") && (given["name"].as_string() == "none"))
        return RecipeValue();
      else if ((kind_of(given) == "name") && (given["name"].as_string() == "auto"))
        automatic = true;
      else
      {
        std::vector<RecipeValue> items;
        if (kind_of(given) == "vector")
          for (size_t k = 0; k < given["elements"].size(); ++k)
            items.push_back(given["elements"][k]);
        else
          items.push_back(given);
        for (const RecipeValue &item : items)
        {
          static const std::string forms = "transform is auto, none, or a list of log(x), logit(x), logit(x, a, b), lower(x, a), upper(x, b)";
          if (kind_of(item) != "call")
            this->error(item, forms + ".");
          std::string f = item["function"].as_string();
          const RecipeValue &arguments = item["arguments"];
          if (((f != "log") && (f != "logit") && (f != "lower") && (f != "upper")) || (arguments.size() < 1) ||
              (kind_of(arguments[0]) != "name"))
            this->error(item, forms + ".");
          std::string x = arguments[0]["name"].as_string();
          if (std::find(moved.begin(), moved.end(), x) == moved.end())
            this->error(arguments[0], "transform: " + x + " is not moved by this move.");
          size_t expected = (f == "log") ? 1 : (f == "logit") ? ((arguments.size() == 1) ? 1 : 3) : 2;
          if (arguments.size() != expected)
            this->error(item, f + ": " + ((f == "log") ? "log(x)" : (f == "logit") ? "logit(x) or logit(x, a, b)" :
                                          (f == "lower") ? "lower(x, a)" : "upper(x, b)") + ".");
          double a = 0.0, b = (f == "logit") ? 1.0 : 0.0;
          if ((f == "logit") && (arguments.size() == 3))
          {
            if (!this->compile_time_number(arguments[1], a) || !this->compile_time_number(arguments[2], b))
              this->error(item, "logit(x, a, b): the bounds are numbers (or constants).");
            if (!(a < b))
              this->error(item, "logit(x, a, b): a must be below b.");
          }
          if ((f == "lower") && !this->compile_time_number(arguments[1], a))
            this->error(item, "lower(x, a): the bound is a number (or a constant).");
          if ((f == "upper") && !this->compile_time_number(arguments[1], b))
            this->error(item, "upper(x, b): the bound is a number (or a constant).");
          list.push_back(transform_item(x, f, a, b));
        }
        return list;
      }
    }
    if (!automatic)
      return RecipeValue();
    for (const std::string &x : moved)
    {
      RecipeValue item = this->automatic_transform(x);
      if (!item.is_null())
        list.push_back(item);
    }
    return (list.size() > 0) ? list : RecipeValue();
  }

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
    this->allow_only(call, arguments, {"proposal", "log", "factors", "cut_factors", "adapt", "smc_adapt", "transform"}, 1);
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
      RecipeValue node = this->builtin_node(proposal, true, "", true);
      variable = node["variables"].as_string();
      entry.push_back(node, p + "_proposal");
      // the scale it moves on (C209): auto by default for hmc
      if (p != "imh")
      {
        RecipeValue transform = this->move_transform(arguments.has("transform") ? &arguments.named["transform"] : NULL,
                                                     variable, type == "hmc");
        if (!transform.is_null())
          entry.push_back(transform, "transform");
      }
      else if (arguments.has("transform"))
        this->error(arguments.named["transform"], kind + ": an independent proposal does not move on a transformed scale.");
      if (arguments.has("log"))
        this->error(arguments.named["log"], kind + ": log goes with a proposal written by the user.");
    }
    else
    {
      if (arguments.has("transform"))
        this->error(arguments.named["transform"], kind + ": transform goes with a built-in proposal "
                    "(a proposal written by the user moves on the scale it chooses).");
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
    else
    {
      // (I4b) a move on known variables evaluates only the factors that may
      // read them: those from model statements reading one of them, and every
      // other factor (what it reads is not known)
      std::vector<size_t> numbers = this->automatic_factor_index(variable);
      if (!numbers.empty())
        entry.push_back(this->placeholder(this->generate_factor_index_numbers(numbers), p + "_factor_index"),
                        p + "_factor_index");
    }
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
  // --- when to stop adapting (until = …, dev/plan.md 3.18 item 10) ------------

  // a number in a criterion (constants by their values when compiling)
  double criterion_number(const RecipeValue &expression, const std::string &what)
  {
    std::map<std::string, double> constants;
    for (const std::string &constant : this->constant_names)
    {
      double value;
      if (this->constant_number(constant, value))
        constants[constant] = value;
    }
    arma::mat number;
    if (!this->numeric(substitute(expression, constants), number) || (number.n_elem != 1))
      this->error(expression, what + " is a number.");
    return number(0);
  }

  // until = criterion, as the recipe's node
  RecipeValue adaptation_criterion(const RecipeValue &expression)
  {
    static const std::string forms = "iterations(n), stable(tolerance, window, watch, change), acceptance(tolerance, window), "
                                     "all(…) or any(…)";
    if (kind_of(expression) != "call")
      this->error(expression, "until is a criterion: " + forms + ".");
    std::string f = expression["function"].as_string();
    Arguments arguments = this->split_arguments(expression);
    RecipeValue node = RecipeValue::list();
    node.push_back(RecipeValue(f), "criterion");
    auto argument = [&](size_t position, const std::string &name) -> RecipeValue {
      if (arguments.has(name))
        return arguments.named[name];
      return (position < arguments.positional.size()) ? arguments.positional[position] : RecipeValue();
    };
    if (f == "iterations")
    {
      this->allow_only(expression, arguments, {"n"}, 1);
      RecipeValue n = argument(0, "n");
      if (n.is_null())
        this->error(expression, "iterations(n): the number of adaptation steps is missing.");
      double value = this->criterion_number(n, "iterations(n): n");
      if ((value < 1) || (value != std::floor(value)))
        this->error(n, "iterations(n): n is a whole number from 1.");
      node.push_back(RecipeValue(value), "n");
      return node;
    }
    if ((f == "stable") || (f == "acceptance"))
    {
      bool stable = (f == "stable");
      if (stable)
        this->allow_only(expression, arguments, {"tolerance", "window", "watch", "change"}, 2);
      else
        this->allow_only(expression, arguments, {"tolerance", "window"}, 2);
      // defaults: placeholders to be set from tests (dev/plan.md 3.18)
      RecipeValue tolerance = argument(0, "tolerance"), window = argument(1, "window");
      double tolerance_value = tolerance.is_null() ? (stable ? 0.01 : 0.05) : this->criterion_number(tolerance, f + ": tolerance");
      double window_value = window.is_null() ? 100.0 : this->criterion_number(window, f + ": window");
      if (!(tolerance_value >= 0.0))
        this->error(tolerance, f + ": the tolerance is at least 0.");
      if ((window_value < 1) || (window_value != std::floor(window_value)))
        this->error(window.is_null() ? expression : window, f + ": the window is a whole number from 1.");
      node.push_back(RecipeValue(tolerance_value), "tolerance");
      node.push_back(RecipeValue(window_value), "window");
      if (stable && arguments.has("watch"))
      {
        static const std::set<std::string> watchable = {"default", "scale", "covariance", "covariance_elements", "mean",
                                                        "component_scales"};
        const RecipeValue &watch = arguments.named["watch"];
        std::vector<std::string> names;
        auto one = [&](const RecipeValue &e) {
          if ((kind_of(e) != "name") || !watchable.count(e["name"].as_string()))
            this->error(e, "stable: watch is default, scale, covariance, covariance_elements, mean or component_scales, "
                        "or a list of these.");
          names.push_back(e["name"].as_string());
        };
        if (kind_of(watch) == "vector")
          for (size_t k = 0; k < watch["elements"].size(); ++k)
            one(watch["elements"][k]);
        else
          one(watch);
        node.push_back(RecipeValue(names), "watch");
      }
      if (stable && arguments.has("change"))
      {
        const RecipeValue &change = arguments.named["change"];
        if ((kind_of(change) != "name") || ((change["name"].as_string() != "relative") && (change["name"].as_string() != "absolute")))
          this->error(change, "stable: change is relative or absolute.");
        node.push_back(RecipeValue(change["name"].as_string()), "change");
      }
      return node;
    }
    if ((f == "all") || (f == "any"))
    {
      this->allow_only(expression, arguments, {}, 1000);
      if (arguments.positional.empty())
        this->error(expression, f + "(…) holds one or more criteria.");
      RecipeValue items = RecipeValue::list();
      for (const RecipeValue &item : arguments.positional)
        items.push_back(this->adaptation_criterion(item));
      node.push_back(items, "items");
      return node;
    }
    this->error(expression, "until is a criterion: " + forms + ".");
  }

  // the adaptor's call without until = (lowered separately into `until`)
  RecipeValue without_until(const RecipeValue &adaptor, RecipeValue &until)
  {
    if ((kind_of(adaptor) != "call") || !adaptor.contains("arguments"))
      return adaptor;
    const RecipeValue &arguments = adaptor["arguments"];
    std::vector<std::string> keys = arguments.names();
    RecipeValue kept = RecipeValue::list();
    for (size_t i = 0; i < arguments.size(); ++i)
    {
      std::string key = (i < keys.size()) ? keys[i] : "";
      if (key == "until")
        until = this->adaptation_criterion(arguments[i]);
      else
        kept.push_back(arguments[i], key);
    }
    if (until.is_null())
      return adaptor;
    // a new node (values share their nodes: the original is not changed)
    RecipeValue call = RecipeValue::list();
    std::vector<std::string> fields = adaptor.names();
    for (size_t i = 0; i < adaptor.size(); ++i)
      call.push_back((fields[i] == "arguments") ? kept : adaptor[i], fields[i]);
    return call;
  }

  void lower_adaptor(const RecipeValue &adaptor_with_until, const std::string &p, const std::string &variable, RecipeValue &entry)
  {
    RecipeValue until;
    RecipeValue adaptor = this->without_until(adaptor_with_until, until);
    if (!until.is_null())
    {
      if (this->builtin_name(adaptor) == "averaged")
        this->error(adaptor_with_until, "until goes with a built-in adaptor (e.g. scaled_am(until = iterations(1000))).");
      entry.push_back(until, p + "_adaptor_until");
    }
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
    RecipeValue mfds_transform;
    RecipeValue algorithm_transform;

    auto method = [&](const std::string &key, const RecipeValue &node) {
      RecipeValue item = RecipeValue::list();
      item.push_back(node, key);
      methods.push_back(item);
    };

    const RecipeValue &statements = section["statements"];
    for (size_t i = 0; i < statements.size(); ++i)
    {
      const RecipeValue &statement = statements[i];
      if (statement.contains("keyword") && statement["keyword"].as_string() == "for")
      {
        // a scan of moves on elements (I4c)
        Move move = this->element_move_loop(statement);
        for (const Move &m : moves)
          if (m.name == move.name)
            this->error(statement, "move " + move.name + " is defined more than once.");
        moves.push_back(move);
        continue;
      }
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
        if ((kind_of(value) == "call") && (value["function"].as_string() == "gibbs"))
          // gibbs(mu) on a scalar, as a move on its one element (C206)
          moves.push_back(this->element_move_loop(statement));
        else
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
      else if ((target == "transform") && !this->is_function_reference(value))
      {
        // the built-in transforms for the whole algorithm (IENKI, C209):
        // auto (every parameter with a support in the model statements), none,
        // or a list [log(x), logit(p, a, b), …] of any parameters
        if (!algorithm_transform.is_null() || !mfds_transform.is_null())
          this->error(statement, "\"transform\" is set more than once.");
        std::string parameters;
        for (const std::string &name : this->model_parameter_names())
          parameters = parameters + (parameters.empty() ? "" : ";") + name;
        algorithm_transform = this->move_transform(&value, parameters, false);
        if (algorithm_transform.is_null())
          algorithm_transform = RecipeValue::list(); // (none, or nothing to transform)
      }
      else if (target == "transform")
      {
        // MFDS: the map from the Gaussian scale the ensemble moves on to the
        // model's parameters, under the same names (dev/plan.md 3.16); the
        // variables it does not give keep their values
        if (value["namespace"].as_string() != "cpp")
          this->error(value, "transform is a C++ function (cpp::f) giving the model's parameters from the Gaussian scale.");
        std::string f = this->cpp_function(value, "transform");
        std::string wrapper = "ilike_mfds_transform_" + std::to_string(++this->generated_count);
        this->generated_code += "Parameters " + wrapper + "(const Parameters &u)\n{\n  Parameters theta = " + f +
                                "(u);\n  return theta.merge(u);\n}\n\n";
        mfds_transform = this->placeholder(wrapper, "mfds_transform");
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
    if (!mfds_transform.is_null())
      algorithm.push_back(mfds_transform, "mfds_transform");
    if (!algorithm_transform.is_null() && (algorithm_transform.size() > 0))
      algorithm.push_back(algorithm_transform, "algorithm_transform");
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

RecipeValue lower_ilike_language(const RecipeValue &trees, const RecipeValue &constant_values)
{
  Lowering lowering(trees);
  lowering.set_constant_values(constant_values);
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
