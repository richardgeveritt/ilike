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
    if (this->unnormalised_entries.size() > 0)
    {
      this->recipe.push_back(this->unnormalised_entries, "unnormalised");
      std::vector<std::string> names;
      for (const auto &kv : this->normaliser_parameters)
        names.push_back(kv.first);
      if (!names.empty())
        this->recipe.push_back(RecipeValue(names), "unnormalised_parameters");
    }

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
    output.push_back(RecipeValue(this->notes), "notes");
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
  // (dev/plan.md 3.21, step (b)) integrate u with is(…): the variable, the
  // method's call, the statement (for messages); its factor is named u
  struct Integration
  {
    std::string variable;
    RecipeValue method;
    RecipeValue node;
  };
  std::vector<Integration> integrations;
  std::set<std::string> integrated_variables;
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
  // gradient = autodiff (dev/plan.md 3.10): the factors whose gradients are
  // finite differences, and the notes compile() prints
  // unnormalised distributions (dev/plan.md 3.20, S3): the statements, and
  // the parameters their normalisers depend on (parameter -> the variable
  // whose distribution it is a parameter of)
  struct UnnormalisedStatement
  {
    std::string variable, kind, factor;
    std::vector<std::string> normaliser;
  };
  std::vector<UnnormalisedStatement> unnormalised;
  std::map<std::string, std::string> normaliser_parameters;
  RecipeValue unnormalised_entries = RecipeValue::list();
  std::vector<std::string> finite_difference_gradients;
  bool finite_difference_helper_emitted = false;
  std::map<std::string, std::string> template_log_wrappers; // "f prior"/"f likelihood" -> double instance

  std::vector<std::string> notes;

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
  // gradient = autodiff, autodiff(finite_difference) or autodiff(forward)
  // (dev/plan.md 3.10, S-5): the backend asked for ("auto" for the first), or
  // "" if the value is not autodiff
  std::string autodiff_backend(const RecipeValue &value) const
  {
    std::string kind = kind_of(value);
    if ((kind == "name") && value["namespace"].as_string().empty() && (value["name"].as_string() == "autodiff"))
      return "auto";
    if ((kind != "call") || !value["namespace"].as_string().empty() || (value["function"].as_string() != "autodiff"))
      return "";
    const RecipeValue &arguments = value["arguments"];
    std::vector<std::string> keys = arguments.names();
    if ((arguments.size() == 1) && (keys.empty() || keys[0].empty()) && (kind_of(arguments[0]) == "name") &&
        arguments[0]["namespace"].as_string().empty())
    {
      std::string backend = arguments[0]["name"].as_string();
      if ((backend == "finite_difference") || (backend == "forward") || (backend == "reverse"))
        return backend;
    }
    this->error(value, "autodiff takes no argument (compile() chooses how), or one of finite_difference, forward, reverse.");
    return ""; // (not reached)
  }

  // central differences of a log density, emitted once into the generated
  // code: each element of the variable moved by h = ε^(1/3) max(1, |x|)
  // (the h balancing the truncation error, O(h²), against the rounding
  // error, O(ε/h)); the parameters are copied deeply, since a copy of
  // Parameters shares the variables' memory
  void emit_finite_difference_helper()
  {
    if (this->finite_difference_helper_emitted)
      return;
    this->finite_difference_helper_emitted = true;
    this->generated_code +=
      "template <class F>\n"
      "static arma::mat ilike_finite_difference_gradient(F f, const std::string &variable, const Parameters &parameters)\n"
      "{\n"
      "  Parameters moved = parameters.deep_copy();\n"
      "  arma::mat &x = moved[variable];\n"
      "  arma::mat gradient(x.n_rows, x.n_cols);\n"
      "  const double scale = std::cbrt(std::numeric_limits<double>::epsilon());\n"
      "  for (arma::uword k = 0; k < x.n_elem; ++k)\n"
      "  {\n"
      "    double value = x(k);\n"
      "    double h = scale * std::max(1.0, std::abs(value));\n"
      "    double up = value + h, down = value - h;\n"
      "    x(k) = up;\n"
      "    double f_up = f(moved);\n"
      "    x(k) = down;\n"
      "    double f_down = f(moved);\n"
      "    x(k) = value;\n"
      "    gradient(k) = (f_up - f_down) / (up - down); // (up - down: the step actually taken)\n"
      "  }\n"
      "  return gradient;\n"
      "}\n\n";
  }

  // the gradient function generated for gradient = autodiff (dev/plan.md
  // 3.10 stage 1): for now, central differences of the log density, always
  // in C++: of the user's function, or, for a log density in another
  // language, of the C++ wrapper the R binding generates for it (placeholder
  // {ilike_finite_difference, log: the log's placeholder, role}); either way
  // an ordinary gradient function, so the engine sees no difference
  // a log density written in C++ as a template on the number type,
  // template <class T> T f(const ParametersT<T> &parameters[, const Data &data])
  // (dev/plan.md 3.10 stage 2; ilike_autodiff.h)
  bool is_template_log(const RecipeValue &log)
  {
    if (this->is_foreign(log) || !this->is_function_reference(log))
      return false;
    std::string f = this->cpp_function(log, "log");
    return find_cpp_function_signature(this->user_code, f).find("ParametersT") != std::string::npos;
  }

  // the double instance of a templated log density, as the ordinary function
  // the log density's role takes; ParametersT<double> is a view on the
  // Parameters, so it costs what the function written with Parameters costs
  std::string template_log_wrapper(const RecipeValue &log, bool prior)
  {
    std::string f = this->cpp_function(log, "log");
    std::string key = f + (prior ? " prior" : " likelihood");
    auto found = this->template_log_wrappers.find(key);
    if (found != this->template_log_wrappers.end())
      return found->second;
    bool takes_data = find_cpp_function_signature(this->user_code, f).find("Data") != std::string::npos;
    std::string w = "ilike_generated_log_" + std::to_string(++this->generated_count);
    this->generated_code += "double " + w + "(const Parameters &parameters" + (prior ? "" : ", const Data &data") +
                            ")\n{\n  return ::" + f + "<double>(ParametersT<double>(parameters)" +
                            (takes_data ? ", data" : "") + ");\n}\n\n";
    this->template_log_wrappers[key] = w;
    return w;
  }

  // the placeholder of a log density: the function, or a templated
  // function's double instance
  RecipeValue log_placeholder(const RecipeValue &log, bool prior, const std::string &what)
  {
    std::string role = std::string("evaluate_log_") + (prior ? "prior" : "likelihood");
    if (this->is_template_log(log))
      return this->placeholder(this->template_log_wrapper(log, prior), role);
    return this->function_placeholder(log, role, what);
  }

  RecipeValue generated_gradient(const RecipeValue &where, const std::string &backend, Arguments &arguments, bool prior,
                                 const std::string &name)
  {
    std::string suffix = prior ? "prior" : "likelihood";
    if (!arguments.has("log"))
      this->error(where, "gradient = autodiff differentiates the log density, so the log density is needed too: "
                         "density(log = …, gradient = autodiff).");
    const RecipeValue &log = arguments.named["log"];
    bool template_log = this->is_template_log(log);
    if (((backend == "forward") || (backend == "reverse")) && !template_log)
      this->error(where, "autodiff(" + backend + ") needs the log density written in C++ as a template on the number type: "
                         "template <class T> T f(const ParametersT<T> &parameters, const Data &data) (see ilike_autodiff.h).");
    if (template_log && (backend != "finite_difference"))
    {
      // exact: forward mode (one Dual evaluation per element) or reverse mode
      // (one Var evaluation recorded on a tape, then one sweep back, whatever
      // the number of elements, but each operation costs more); by default
      // forward for fewer than 10 elements, reverse otherwise (C239: the two
      // took the same time at 10 elements on a logistic regression, forward
      // 3 times faster at 2, reverse 3.5 times faster at 50)
      std::string f = this->cpp_function(log, "log");
      bool takes_data = find_cpp_function_signature(this->user_code, f).find("Data") != std::string::npos;
      auto call = [&](const std::string &number) {
        return std::string(number == "Dual" ? "forward_gradient" : "reverse_gradient") + "([&](const ParametersT<" + number +
               "> &p) { return ::" + f + "<" + number + ">(p" + (takes_data ? ", data" : "") + "); }, variable, parameters)";
      };
      std::string body;
      if (backend == "forward")
        body = "  return " + call("Dual") + ";\n";
      else if (backend == "reverse")
        body = "  return " + call("Var") + ";\n";
      else
        body = "  if (parameters[variable].n_elem < 10)\n    return " + call("Dual") + ";\n  return " + call("Var") + ";\n";
      std::string g = "ilike_generated_gradient_" + std::to_string(++this->generated_count);
      this->generated_code += "arma::mat " + g + "(const std::string &variable, const Parameters &parameters" +
                              (prior ? "" : ", const Data &data") + ")\n{\n" + body + "}\n\n";
      return this->placeholder(g, "evaluate_gradient_log_" + suffix);
    }
    this->finite_difference_gradients.push_back(name);
    this->notes.push_back(std::string("The gradient of the ") + (prior ? "prior " : "factor ") + name +
                          " is computed by finite differences of its log density, so it is approximate: "
                          "metropolis_hastings moves using it stay exact, unadjusted moves refuse it" +
                          std::string(template_log || this->is_foreign(log) ? "." :
                                      " (a log density written in C++ as a template on the number type gets an exact gradient)."));
    this->emit_finite_difference_helper();
    if (this->is_foreign(log))
    {
      RecipeValue output = RecipeValue::list();
      output.push_back(RecipeValue(prior ? "prior" : "likelihood"), "ilike_finite_difference");
      output.push_back(this->foreign_placeholder(log, "evaluate_log_" + suffix, ""), "log");
      output.push_back(RecipeValue("evaluate_gradient_log_" + suffix), "role");
      return output;
    }
    // (a templated log density is differenced through its double instance)
    std::string f = template_log ? this->template_log_wrapper(log, prior) : this->cpp_function(log, "log");
    bool takes_data = template_log ? !prior : (find_cpp_function_signature(this->user_code, f).find("Data") != std::string::npos);
    std::string g = "ilike_generated_gradient_" + std::to_string(++this->generated_count);
    // (::f: the generated function's arguments must not hide the user's function)
    this->generated_code += "arma::mat " + g + "(const std::string &variable, const Parameters &parameters" +
                            (prior ? "" : ", const Data &data") + ")\n{\n" +
                            "  return ilike_finite_difference_gradient([&](const Parameters &p) { return (double)::" + f +
                            "(p" + (takes_data ? ", data" : "") + "); }, variable, parameters);\n}\n\n";
    return this->placeholder(g, "evaluate_gradient_log_" + suffix);
  }

  RecipeValue density_functions(const RecipeValue &value, bool prior, const std::string &name)
  {
    RecipeValue factor = RecipeValue::list();
    std::string suffix = prior ? "prior" : "likelihood";
    if (this->is_function_value(value))
    {
      factor.push_back(this->log_placeholder(value, prior, "the log density"), "evaluate_log_" + suffix);
      return factor;
    }
    if (kind_of(value) != "call" || value["function"].as_string() != "density")
      this->error(value, std::string("a ") + (prior ? "prior" : "factor") +
                  " is a function (cpp::f) or density(log = cpp::f, …).");
    Arguments arguments = this->split_arguments(value);
    std::set<std::string> allowed = {"log", "gradient", "second_gradient"};
    if (prior)
      allowed.insert("simulate");
    else
      allowed.insert("maximum"); // (dev/plan.md 3.21)
    this->allow_only(value, arguments, allowed, 0);
    if (arguments.has("log"))
      factor.push_back(this->log_placeholder(arguments.named["log"], prior, "log"), "evaluate_log_" + suffix);
    if (arguments.has("simulate"))
      factor.push_back(this->function_placeholder(arguments.named["simulate"], "simulate_prior", "simulate"),
                       "simulate_prior");
    // maximum = cpp::f: the log likelihood's largest value over the state, at
    // the data (a function with the log likelihood's arguments), so that a
    // particle filter can stop early (dev/plan.md 3.21)
    if (arguments.has("maximum"))
      factor.push_back(this->function_placeholder(arguments.named["maximum"], "evaluate_log_likelihood", "maximum"),
                       "maximum_log_likelihood");
    if (arguments.has("gradient"))
    {
      std::string backend = this->autodiff_backend(arguments.named["gradient"]);
      if (backend.empty())
        factor.push_back(this->function_placeholder(arguments.named["gradient"], "evaluate_gradient_log_" + suffix, "gradient"), "evaluate_gradient_log_" + suffix);
      else
        factor.push_back(this->generated_gradient(arguments.named["gradient"], backend, arguments, prior, name),
                         "evaluate_gradient_log_" + suffix);
    }
    if (arguments.has("second_gradient") && !this->autodiff_backend(arguments.named["second_gradient"]).empty())
      this->error(arguments.named["second_gradient"], "second_gradient = autodiff is not written yet.");
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

      if (keyword == "integrate")
      {
        // integrate u with is(points = N) (dev/plan.md 3.21, step (b)): the
        // factor of the data that depend on u, estimated by an importance
        // sampler over u for each value of the other parameters; the factor
        // is named after u and sits here (finish_model_statements fills it)
        if (op != "with" || kind_of(value) != "call" || !value["namespace"].as_string().empty() ||
            (value["function"].as_string() != "is" && value["function"].as_string() != "abc"))
          this->error(statement, "integrate " + name + " with is(points = N) (or abc(points = N)): a variable "
                      "integrated out by importance sampling.");
        if (!this->integrations.empty())
          this->error(statement, "only one variable can be integrated out, for now (it may be a vector or matrix).");
        this->add_factor(statement, name, RecipeValue::list());
        this->integrations.push_back({name, value, statement});
        this->integrated_variables.insert(name);
        this->mfds_reasons.push_back(name + " is integrated out by an inner algorithm");
        continue;
      }

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
        this->add_factor(statement, name, this->density_functions(value, true, name));
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
          this->add_factor(statement, name, this->density_functions(value, false, name));
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
          if (keyword == "for" || op == "<-" || op == "+=" || keyword == "integrate")
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

  // an expression as ilike text (for statements the lowering writes itself)
  std::string ilike_text(const RecipeValue &expression) const
  {
    std::string kind = kind_of(expression);
    if (kind == "number")
      return number_text(expression["value"].as_double());
    if (kind == "name")
      return (expression["namespace"].as_string().empty() ? std::string("") : expression["namespace"].as_string() + "::") +
             expression["name"].as_string();
    if (kind == "index")
    {
      const RecipeValue &indices = expression["indices"];
      std::string text = this->ilike_text(expression["target"]) + "[";
      for (size_t k = 0; k < indices.size(); ++k)
        text += (k ? ", " : "") + this->ilike_text(indices[k]);
      return text + "]";
    }
    if (kind == "unary")
      return "(" + expression["operator"].as_string() + this->ilike_text(expression["operand"]) + ")";
    if (kind == "binary")
      return "(" + this->ilike_text(expression["left"]) + " " + expression["operator"].as_string() + " " +
             this->ilike_text(expression["right"]) + ")";
    if (kind == "call")
    {
      const RecipeValue &arguments = expression["arguments"];
      std::vector<std::string> keys = arguments.names();
      std::string text = this->ilike_text_name(expression) + "(";
      for (size_t k = 0; k < arguments.size(); ++k)
        text += (k ? ", " : "") + ((k < keys.size() && !keys[k].empty()) ? keys[k] + " = " : std::string("")) +
                this->ilike_text(arguments[k]);
      return text + ")";
    }
    this->error(expression, "this argument cannot be used here (a number, a name, an element or an arithmetic expression).");
    return ""; // (not reached)
  }
  std::string ilike_text_name(const RecipeValue &call) const
  {
    return (call["namespace"].as_string().empty() ? std::string("") : call["namespace"].as_string() + "::") +
           call["function"].as_string();
  }

  static bool is_unnormalised_builtin(const RecipeValue &value)
  {
    if (kind_of(value) != "call")
      return false;
    std::string name_space = value["namespace"].as_string();
    std::string f = value["function"].as_string();
    return (name_space.empty() || (name_space == "ilike")) && ((f == "ising") || (f == "potts"));
  }

  // x ~ ising(beta, nb, nn) and x ~ potts(beta, q, nb, nn) (dev/plan.md
  // 3.20, S3): written as their base measure (x[i] ~ bern(0.5), or
  // categorical with equal probabilities) and their potential, a factor of
  // terms named x_potential, so that moves on the states recompute only the
  // neighbouring terms, the prior's simulation draws the base measure, and
  // annealing the potential starts from it. The normaliser Z(beta) is not
  // part of either: the parameters it depends on are recorded, and the moves
  // and algorithms that would need Z are refused. Also generated: the
  // simulators of x given beta (heat-bath sweeps; coupling from the past for
  // ising).
  void expand_unnormalised(const RecipeValue &statement, const std::vector<ModelLoop> &loops)
  {
    const RecipeValue &value = statement["value"];
    std::string kind = value["function"].as_string();
    if (!loops.empty())
      this->error(statement, kind + "(…) is the distribution of a whole field: write it outside loops (x ~ " + kind + "(…)).");
    if (kind_of(statement["target"]) != "name")
      this->error(statement, kind + "(…)'s target is a whole variable (x ~ " + kind + "(…), without an index).");
    const std::string x = statement["target"]["name"].as_string();
    Arguments arguments = this->split_arguments(value);
    std::vector<std::string> order = (kind == "ising") ? std::vector<std::string>{"beta", "nb", "nn"}
                                                       : std::vector<std::string>{"beta", "q", "nb", "nn"};
    std::set<std::string> allowed(order.begin(), order.end());
    allowed.insert("n");
    allowed.insert("normaliser");
    this->allow_only(value, arguments, allowed, order.size());
    std::map<std::string, RecipeValue> given;
    for (size_t k = 0; k < arguments.positional.size(); ++k)
      given[order[k]] = arguments.positional[k];
    for (const auto &kv : arguments.named)
      given[kv.first] = kv.second;
    RecipeValue normaliser_method;
    if (given.count("normaliser"))
    {
      normaliser_method = given["normaliser"];
      given.erase("normaliser");
    }
    for (const std::string &name : order)
      if (!given.count(name))
        this->error(value, kind + "(…): the argument " + name + " is missing (" + kind + "(" + joined(order, ", ") + ")).");
    for (const std::string &name : {std::string("nb"), std::string("nn")})
      if ((kind_of(given[name]) != "name") || !this->data_names.count(given[name]["name"].as_string()))
        this->error(given[name], kind + "(…): " + name + " is a data variable (the neighbours' lists, as in sum(k in 1:nn[i], …)).");
    const std::string nb = given["nb"]["name"].as_string(), nn = given["nn"]["name"].as_string();
    // the number of sites: n, or the length of nn given as values
    std::string sites;
    double sites_value = -1.0;
    if (given.count("n"))
      sites = this->ilike_text(given["n"]);
    else
    {
      for (const auto &item : this->literal_data)
        if (item.first == nn)
          sites_value = double(item.second.as_matrix().n_elem);
      if (sites_value < 0.0)
        this->error(value, kind + "(…): the number of sites is the length of " + nn + " when it is given as values in the data section; "
                    "otherwise give it, " + kind + "(…, n = …).");
      sites = number_text(sites_value);
    }
    double q = 2.0;
    if (kind == "potts")
    {
      arma::mat number;
      if (!this->numeric(given["q"], number) || (number.n_elem != 1) || (number(0) < 2.0) || (number(0) != std::floor(number(0))))
        this->error(given["q"], "potts(…): q, the number of colours, is a whole number of at least 2, known when compiling.");
      q = number(0);
    }
    std::string beta = this->ilike_text(given["beta"]);
    std::string factor = x + "_potential";

    // the statements, written in ilike and parsed
    std::string base, potential;
    if (kind == "ising")
    {
      base = "bern(0.5)";
      potential = "(" + beta + ") / 2 * (2 * " + x + "[ilike_site] - 1) * sum(ilike_k in 1:" + nn + "[ilike_site], 2 * " + x + "[" + nb +
                  "[ilike_site, ilike_k]] - 1)";
    }
    else
    {
      base = "categorical([";
      for (int c = 0; c < int(q); ++c)
        base += (c ? ", " : "") + number_text(1.0 / q);
      base += "])";
      potential = "(" + beta + ") / 2 * sum(ilike_k in 1:" + nn + "[ilike_site], equals(" + x + "[ilike_site], " + x + "[" + nb +
                  "[ilike_site, ilike_k]]))";
    }
    std::string text = "ilike 1.0\nmodel {\n  for (ilike_site in 1:" + sites + ") " + x + "[ilike_site] ~ " + base + "\n" +
                       "  for (ilike_site in 1:" + sites + ") factor " + factor + " += " + potential + "\n}\n";
    RecipeValue tree = parse_ilike_language(text, this->file + " (" + kind + " of " + x + ")");
    const RecipeValue &statements = tree["sections"][0]["statements"];
    for (size_t k = 0; k < statements.size(); ++k)
      this->collect_model_statement(statements[k], std::vector<ModelLoop>());

    // the parameters the normaliser depends on
    UnnormalisedStatement record;
    record.variable = x;
    record.kind = kind;
    record.factor = factor;
    std::set<std::string> found;
    this->read_parameters(given["beta"], found);
    if (given.count("n"))
      this->read_parameters(given["n"], found);
    for (const std::string &parameter : found)
    {
      record.normaliser.push_back(parameter);
      if (normaliser_method.is_null())
        this->normaliser_parameters[parameter] = x;
    }
    this->unnormalised.push_back(record);

    // the simulators of x given the rest
    CodeContext context;
    std::string beta_code = this->scalar_code(given["beta"], context);
    std::string sites_code = given.count("n") ? this->scalar_code(given["n"], context) : number_text(sites_value);
    std::string sweep = "ilike_unnormalised_sweep_" + x, exact = "ilike_unnormalised_exact_" + x;
    std::string head =
      "  const arma::mat &nb = data[" + quoted(nb) + "];\n"
      "  const arma::mat &nn = data[" + quoted(nn) + "];\n"
      "  const double beta = " + beta_code + ";\n"
      "  const size_t n = (size_t)std::lround(" + sites_code + ");\n";
    std::string code = "// " + kind + " field " + x + " (dev/plan.md 3.20, S3): sweeps of single-site heat-bath updates given beta\n"
                       "void " + sweep + "(RandomNumberGenerator &rng, Parameters &parameters, const Data &data)\n{\n" + head +
                       "  arma::mat &x = parameters[" + quoted(x) + "];\n"
                       "  if (x.n_elem != n)\n  {\n    // (a first draw: from the base measure)\n    x = arma::mat(n, 1);\n";
    if (kind == "ising")
      code += "    for (size_t i = 0; i < n; ++i) x(i) = (runif(rng) < 0.5) ? 1.0 : 0.0;\n  }\n"
              "  for (size_t i = 0; i < n; ++i)\n  {\n"
              "    double s = 0.0;\n"
              "    for (size_t k = 0; k < (size_t)nn(i); ++k) s += 2.0 * x((size_t)nb(i, k) - 1) - 1.0;\n"
              "    x(i) = (runif(rng) < 1.0 / (1.0 + std::exp(-2.0 * beta * s))) ? 1.0 : 0.0;\n  }\n}\n\n";
    else
      code += "    for (size_t i = 0; i < n; ++i) x(i) = 1.0 + std::min(std::floor(runif(rng) * (double)" + number_text(q) + "), (double)" +
              number_text(q - 1.0) + ");\n  }\n"
              "  std::vector<double> w(" + number_text(q) + ");\n"
              "  for (size_t i = 0; i < n; ++i)\n  {\n"
              "    std::fill(w.begin(), w.end(), 0.0);\n"
              "    for (size_t k = 0; k < (size_t)nn(i); ++k) w[(size_t)x((size_t)nb(i, k) - 1) - 1] += 1.0;\n"
              "    double top = *std::max_element(w.begin(), w.end()), total = 0.0;\n"
              "    for (double &v : w) { v = std::exp(beta * (v - top)); total += v; }\n"
              "    double u = runif(rng) * total;\n"
              "    size_t c = 0;\n"
              "    while ((c + 1 < w.size()) && (u > w[c])) { u -= w[c]; ++c; }\n"
              "    x(i) = double(c + 1);\n  }\n}\n\n";
    RecipeValue entry = RecipeValue::list();
    entry.push_back(RecipeValue(x), "variable");
    entry.push_back(RecipeValue(kind), "kind");
    entry.push_back(this->placeholder(sweep, "unnormalised_sweep"), "sweep");
    if (kind == "ising")
    {
      // monotone coupling from the past (Propp & Wilson 1996): two chains
      // started from all 0 and all 1 at time -T, run with the same random
      // numbers; heat-bath updates are monotone for beta >= 0, so when they
      // meet at time 0, every chain started at -T would have, and the state
      // is an exact draw. T doubles, reusing the random numbers of the most
      // recent sweeps.
      code += "// ising field " + x + ": an exact draw by monotone coupling from the past\n"
              "Parameters " + exact + "(RandomNumberGenerator &rng, const Parameters &parameters, const Data &data)\n{\n" + head +
              "  // (beta < 0: not monotone, so no draw; a field of NaN, which the\n"
              "  // estimate of the normalising constant takes as density 0)\n"
              "  if (beta < 0.0)\n"
              "    return Parameters(" + quoted(x) + ", arma::mat(n, 1, arma::fill::value(arma::datum::nan)));\n"
              "  std::vector<std::vector<double>> uniforms; // (the most recent sweep first: earlier sweeps are appended)\n"
              "  size_t T = 1;\n"
              "  arma::mat upper(n, 1), lower(n, 1);\n"
              "  while (true)\n  {\n"
              "    while (uniforms.size() < T)\n    {\n"
              "      std::vector<double> u(n);\n"
              "      for (size_t i = 0; i < n; ++i) u[i] = runif(rng);\n"
              "      uniforms.push_back(u);\n    }\n"
              "    upper.fill(1.0);\n    lower.fill(0.0);\n"
              "    for (size_t back = T; back > 0; --back)\n"
              "    {\n"
              "      const std::vector<double> &u = uniforms[back - 1]; // (from time -T to time 0)\n"
              "      for (size_t i = 0; i < n; ++i)\n      {\n"
              "        double su = 0.0, sl = 0.0;\n"
              "        for (size_t k = 0; k < (size_t)nn(i); ++k)\n        {\n"
              "          su += 2.0 * upper((size_t)nb(i, k) - 1) - 1.0;\n"
              "          sl += 2.0 * lower((size_t)nb(i, k) - 1) - 1.0;\n        }\n"
              "        upper(i) = (u[i] < 1.0 / (1.0 + std::exp(-2.0 * beta * su))) ? 1.0 : 0.0;\n"
              "        lower(i) = (u[i] < 1.0 / (1.0 + std::exp(-2.0 * beta * sl))) ? 1.0 : 0.0;\n      }\n"
              "    }\n"
              "    if (arma::approx_equal(upper, lower, \"absdiff\", 0.0))\n"
              "      return Parameters(" + quoted(x) + ", upper);\n"
              "    T = 2 * T;\n"
              "    if (double(T) * double(n) > 1e8)\n"
              "      Rcpp::stop(\"exact simulation of ising(...): coupling from the past did not coalesce: beta is far above the critical value, where exact draws take very long; use simulate = sweeps(m), approximate but in a fixed time.\");\n"
              "  }\n}\n\n";
      entry.push_back(this->placeholder(exact, "unnormalised_exact"), "exact");
    }
    std::vector<std::string> normaliser(record.normaliser.begin(), record.normaliser.end());
    entry.push_back(RecipeValue(normaliser), "normaliser");
    this->unnormalised_entries.push_back(entry);

    // log gamma(field | beta) up to a constant (the potential), for the
    // estimates of the normalising constant (dev/plan.md 3.20, S4)
    std::string potential_function = "ilike_unnormalised_potential_" + x;
    code += "double " + potential_function + "(const arma::mat &field, const Parameters &parameters, const Data &data)\n{\n" + head +
            "  double total = 0.0;\n"
            "  for (size_t i = 0; i < n; ++i)\n    for (size_t k = 0; k < (size_t)nn(i); ++k)\n" +
            ((kind == "ising") ? std::string("      total += (2.0 * field(i) - 1.0) * (2.0 * field((size_t)nb(i, k) - 1) - 1.0);\n")
                               : std::string("      total += (field(i) == field((size_t)nb(i, k) - 1)) ? 1.0 : 0.0;\n")) +
            "  return 0.5 * beta * total;\n}\n\n";
    this->deferred_function_code += code;

    if (!normaliser_method.is_null())
      this->normaliser_factor(statement, normaliser_method, x, kind, record, sweep, exact, potential_function);
  }

  // normaliser = sav(points = M, at = theta_hat, simulate = exact | sweeps(m))
  // (dev/plan.md 3.20, S4): the factor x_normaliser, an estimate of
  // 1/Z(theta) from M draws of the field (UnnormalisedNormaliserEstimator).
  // It reads only the normaliser's parameters, so moves on the field do not
  // evaluate it, and moves on its parameters draw new fields (pseudo-marginal).
  void normaliser_factor(const RecipeValue &statement, const RecipeValue &method, const std::string &x, const std::string &kind,
                         const UnnormalisedStatement &record, const std::string &sweep, const std::string &exact,
                         const std::string &potential_function)
  {
    if ((kind_of(method) != "call") || (method["function"].as_string() != "sav"))
      this->error(method, "normaliser = sav(points = M, at = value, simulate = exact | sweeps(m)): the single (M = 1) or "
                          "multiple auxiliary variable method.");
    if (record.normaliser.empty())
      this->error(method, "the normalising constant of " + x + "'s distribution depends on no parameter, so it cancels: "
                          "normaliser = … is not needed.");
    if (record.normaliser.size() != 1)
      this->error(method, "normaliser = sav(…): the normalising constant of " + x + "'s distribution depends on several parameters (" +
                          joined(record.normaliser, ", ") + "); sav takes one for now.");
    Arguments arguments = this->split_arguments(method);
    this->allow_only(method, arguments, {"points", "at", "simulate"}, 2);
    std::map<std::string, RecipeValue> given;
    std::vector<std::string> order = {"points", "at"};
    for (size_t k = 0; k < arguments.positional.size(); ++k)
      given[order[k]] = arguments.positional[k];
    for (const auto &kv : arguments.named)
      given[kv.first] = kv.second;
    if (!given.count("at"))
      this->error(method, "sav(…): at = value, the value of " + record.normaliser[0] + " at which the auxiliary distribution of the "
                          "draws is taken (a value near the posterior's, e.g. a pseudo-likelihood estimate), is needed.");
    arma::mat number;
    if (!this->numeric(given["at"], number) || (number.n_elem != 1))
      this->error(given["at"], "sav(…): at is a number known when compiling.");
    double at = number(0);
    double points = 1.0;
    if (given.count("points"))
    {
      if (!this->numeric(given["points"], number) || (number.n_elem != 1) || (number(0) < 1.0) || (number(0) != std::floor(number(0))))
        this->error(given["points"], "sav(…): points is a whole number of at least 1.");
      points = number(0);
    }
    // the draws: exact (coupling from the past, ising) by default; sweeps(m)
    // otherwise (approximate: the estimate is then biased, by an amount that
    // falls as m grows)
    bool use_exact = (kind == "ising");
    double sweeps = 100.0;
    if (given.count("simulate"))
    {
      const RecipeValue &simulate = given["simulate"];
      if ((kind_of(simulate) == "name") && (simulate["name"].as_string() == "exact"))
      {
        if (kind != "ising")
          this->error(simulate, "simulate = exact: exact simulation is available for ising(…) fields so far; use sweeps(m).");
        use_exact = true;
      }
      else if ((kind_of(simulate) == "call") && (simulate["function"].as_string() == "sweeps") && (simulate["arguments"].size() == 1) &&
               this->numeric(simulate["arguments"][0], number) && (number.n_elem == 1) && (number(0) >= 1.0))
      {
        use_exact = false;
        sweeps = std::floor(number(0));
      }
      else
        this->error(simulate, "simulate = exact, or simulate = sweeps(m) (m heat-bath sweeps from the base measure, approximate).");
    }
    if (!use_exact)
      this->notes.push_back("The draws of " + x + " for its normalising constant are made by " + number_text(sweeps) +
                            " heat-bath sweeps, so they are approximate, and so is the posterior of " + record.normaliser[0] +
                            " (more sweeps make it closer).");
    RecipeValue settings = RecipeValue::list();
    settings.push_back(RecipeValue(x), "variable");
    settings.push_back(this->placeholder(potential_function, "unnormalised_potential"), "potential");
    settings.push_back(this->placeholder(sweep, "unnormalised_sweep"), "sweep");
    if (use_exact)
      settings.push_back(this->placeholder(exact, "unnormalised_exact"), "exact");
    settings.push_back(RecipeValue(sweeps), "sweeps");
    settings.push_back(RecipeValue(points), "points");
    settings.push_back(RecipeValue(record.normaliser), "normaliser");
    settings.push_back(RecipeValue(std::vector<double>(1, at)), "at");
    RecipeValue factor = RecipeValue::list();
    factor.push_back(settings, "unnormalised_normaliser");
    std::string name = x + "_normaliser";
    this->add_factor(statement, name, factor);
    this->generated_factor_reads[this->factor_numbers.at(name)] =
      std::set<std::string>(record.normaliser.begin(), record.normaliser.end());
  }

  void collect_model_statement(const RecipeValue &statement, std::vector<ModelLoop> loops)
  {
    std::string keyword = statement["keyword"].as_string();
    if ((keyword.empty()) && statement.contains("operator") && (statement["operator"].as_string() == "~") &&
        is_unnormalised_builtin(statement["value"]))
    {
      this->expand_unnormalised(statement, loops);
      return;
    }
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

  // the lookups by name of one statement's code, made before its loops
  struct Hoist
  {
    std::map<std::string, std::string> names; // expression -> local name
    std::string declarations;
  };
  size_t hoist_count = 0;

  struct CodeContext
  {
    bool prior = false;                 // in the prior's functions (no runtime data)
    std::set<std::string> indices;      // loop indices in scope
    std::string data = "data";          // how the runtime data are reached ("data" or "(*data)")
    bool node = false;                  // in a deterministic node's function (data is a pointer, nullptr in the prior)
    Hoist *hoist = nullptr;             // when set, lookups by name are made once, before the statement's loops (hoisted)
  };


  bool is_deterministic(const std::string &name) const
  {
    for (const ModelStatement &st : this->model_statements)
      if (st.op == "<-" && st.variable == name)
        return true;
    return false;
  }

  // a deterministic node computed by a C++ function, x <- cpp::f(a, b, …)
  // (dev/plan.md 3.21 (2)): a whole matrix
  bool is_cpp_node(const std::string &name) const
  {
    if (!this->is_deterministic(name))
      return false;
    const RecipeValue &value = this->deterministic_statement(name).value;
    return (kind_of(value) == "call") && (value["namespace"].as_string() == "cpp");
  }
  std::string cpp_node_call(const std::string &name, const CodeContext &context) const
  {
    std::string data = context.prior ? std::string("nullptr") : (context.node ? std::string("data") : std::string("&data"));
    return "ilike_model_node_matrix_" + name + "(parameters, " + data + ")";
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
      return this->hoisted(context.data + "[" + quoted(name) + "]", "const arma::mat &", context);
    }
    if (this->is_cpp_node(name))
      return this->hoisted(this->cpp_node_call(name, context), "const arma::mat &", context); // (computed once per statement when hoisted)
    if (this->is_deterministic(name))
      this->error(node, name + " is a deterministic node: use it with its indices.");
    return this->hoisted("parameters[" + quoted(name) + "]", "const arma::mat &", context);
  }

  // a lookup by name (data["X"], parameters["beta"], ilike_G.find("beta")),
  // done once before a statement's loops when the context collects them: a
  // lookup builds a string and searches a map, which in the innermost loop
  // (a sum over a row of a matrix, say) cost more than the arithmetic. The
  // declaration (type name = expression) goes to the context's Hoist, which
  // the function's code puts before the statement's loops.
  std::string hoisted(const std::string &expression, const std::string &type, const CodeContext &context)
  {
    if (context.hoist == nullptr)
      return expression;
    auto found = context.hoist->names.find(expression);
    if (found != context.hoist->names.end())
      return found->second;
    std::string name = "ilike_ref" + std::to_string(++this->hoist_count);
    context.hoist->declarations += type + " " + name + " = " + expression + ";\n";
    context.hoist->names[expression] = name;
    return name;
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
    if (kind == "sum")
    {
      CodeContext inner;
      std::string loop = this->sum_loop(expression, context, inner);
      std::string total = "ilike_sum" + std::to_string(++this->sum_temporaries);
      return "([&]() { double " + total + " = 0.0; " + loop + " " + total + " += " +
             this->scalar_code(expression["body"], inner) + "; return " + total + "; }())";
    }
    if (kind == "number")
    {
      // a double literal (C231): 1/2 in a statement is 0.5, not C++'s
      // integer division (whole numbers were written as integer literals)
      std::string text = number_text(expression["value"].as_double());
      if (text.find_first_of(".eEin") == std::string::npos)
        text += ".0";
      return "(" + text + ")";
    }
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
      // (C234) indicators, as BUGS's: equals(a, b) is 1 if a = b, else 0;
      // step(x) is 1 if x >= 0, else 0
      if (f == "equals")
      {
        if (arguments.size() != 2) this->error(expression, "equals takes two arguments.");
        return "(((double)" + argument(0) + " == (double)" + argument(1) + ") ? 1.0 : 0.0)";
      }
      if (f == "step")
      {
        if (arguments.size() != 1) this->error(expression, "step takes one argument.");
        return "(((double)" + argument(0) + " >= 0.0) ? 1.0 : 0.0)";
      }
      this->error(expression, "unknown function " + f + " in a model statement (functions: exp, log, sqrt, abs, pow, "
                  "inv_logit, logit, lgamma, log1p, expm1, sin, cos, tanh, min, max, equals, step).");
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
    if (kind == "sum")
    {
      // (C233) its index is local to it
      this->read_parameters(expression["range"], found, visited);
      std::set<std::string> inner;
      this->read_parameters(expression["body"], inner, visited);
      inner.erase(expression["index"].as_string());
      found.insert(inner.begin(), inner.end());
      return;
    }
    if (!expression.is_list())
      return;
    for (size_t k = 0; k < expression.size(); ++k)
      if (expression[k].is_list())
        this->read_parameters(expression[k], found, visited);
  }

  // (C233) the C++ loop header of sum(k in a:b, ...), and the context of its body
  std::string sum_loop(const RecipeValue &expression, const CodeContext &context, CodeContext &inner)
  {
    const RecipeValue &range = expression["range"];
    if ((kind_of(range) != "binary") || (range["operator"].as_string() != ":"))
      this->error(expression, "sum's range is a:b, as in sum(k in 1:n, x[k]).");
    std::string k = expression["index"].as_string();
    if (context.indices.count(k) || this->data_names.count(k) || this->constant_names.count(k) || this->is_deterministic(k))
      this->error(expression, "sum's index " + k + " is already the name of a loop index, data, a constant or a node.");
    std::set<std::string> found;
    this->read_parameters(range, found);
    for (const std::string &index : context.indices) found.erase(index);
    if (!found.empty())
      this->error(expression, "sum's range depends on the parameter " + *found.begin() + " (it must be known from the data and loops).");
    inner = context;
    inner.indices.insert(k);
    return "for (int " + k + " = (int)std::lround(" + this->scalar_code(range["left"], context) + "); " + k +
           " <= (int)std::lround(" + this->scalar_code(range["right"], context) + "); ++" + k + ")";
  }

  std::string deterministic_call(const std::string &name, const RecipeValue &indices, const RecipeValue &node, const CodeContext &context)
  {
    const ModelStatement &st = this->deterministic_statement(name);
    if (this->is_cpp_node(name))
    {
      std::string matrix = this->variable_matrix(name, node, context);
      if (indices.size() == 0)
        return "(" + matrix + ")(0)";
      std::string checked = "false, " + quoted(this->element_description(name, node));
      if (indices.size() == 1)
        return "ilike_model_at(" + matrix + ", " + this->index_code(indices[0], context) + ", " + checked + ")";
      if (indices.size() == 2)
        return "ilike_model_at(" + matrix + ", " + this->index_code(indices[0], context) + ", " + this->index_code(indices[1], context) + ", " + checked + ")";
      this->error(node, "at most two indices.");
    }
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
      if (!k.empty() && (f != "abc_uniform") && (f != "abc_gaussian")) // (the ABC kernels' optional scale and p are named)
        this->error(call, "model statements take the distribution's arguments in order (no names), as in BUGS.");
    static const std::map<std::string, std::pair<size_t, std::string>> densities = {
      {"norm", {2, "dnorm"}}, {"lnorm", {2, "dlnorm"}}, {"gamma", {2, "dgamma"}}, {"unif", {2, "dunif"}},
      {"exp", {1, "dexp"}}, {"pois", {1, "dpois"}}, {"tnorm", {4, "dtnorm"}}, {"beta", {2, "ilike_model_dbeta"}},
      {"bern", {1, "ilike_model_dbern"}}, {"binom", {2, "ilike_model_dbinom"}},
      {"invgamma", {2, "ilike_model_dinvgamma"}}, {"t", {3, "ilike_model_dt"}}, {"cauchy", {2, "ilike_model_dcauchy"}},
      {"halfnorm", {1, "ilike_model_dhalfnorm"}}, {"halfcauchy", {1, "ilike_model_dhalfcauchy"}}, {"halft", {2, "ilike_model_dhalft"}},
      {"negbinom", {2, "ilike_model_dnbinom"}}, {"geom", {1, "ilike_model_dgeom"}}, {"vonmises", {2, "ilike_model_dvonmises"}},
      {"chisq", {1, "ilike_model_dchisq"}}, {"weibull", {2, "ilike_model_dweibull"}}, {"laplace", {2, "ilike_model_dlaplace"}},
      {"logis", {2, "ilike_model_dlogis"}}};
    // (C227) a vector argument or a whole target
    auto check = [&](size_t n, const std::string &form) {
      if (arguments.size() != n)
        this->error(call, f + " takes " + form + ".");
    };
    std::string vx = "arma::colvec(arma::vectorise(" + x + "))";
    if (f == "categorical")
    {
      check(1, "a vector of probabilities");
      return "ilike_model_dcategorical(" + x + ", arma::colvec(arma::vectorise(" + this->matrix_code(arguments[0], context) + ")))";
    }
    if (f == "multinom")
    {
      check(2, "the number of trials and a vector of probabilities");
      return "ilike_model_dmultinom(" + vx + ", " + this->scalar_code(arguments[0], context) + ", arma::colvec(arma::vectorise(" +
             this->matrix_code(arguments[1], context) + ")))";
    }
    if (f == "dirichlet")
    {
      check(1, "a vector of concentrations");
      return "ilike_model_ddirichlet(" + vx + ", arma::colvec(arma::vectorise(" + this->matrix_code(arguments[0], context) + ")))";
    }
    if (f == "mvt")
    {
      check(3, "the degrees of freedom, a location vector and a scale matrix");
      return "ilike_model_dmvt(" + vx + ", " + this->scalar_code(arguments[0], context) + ", arma::colvec(arma::vectorise(" +
             this->matrix_code(arguments[1], context) + ")), " + this->matrix_code(arguments[2], context) + ")";
    }
    if ((f == "wishart") || (f == "invwishart"))
    {
      check(2, "the degrees of freedom and a scale matrix");
      return "ilike_model_d" + f + "(arma::mat(" + x + "), " + this->scalar_code(arguments[0], context) + ", " +
             this->matrix_code(arguments[1], context) + ")";
    }
    if (f == "lkj")
    {
      if ((arguments.size() != 1) && (arguments.size() != 2))
        this->error(call, "lkj takes the shape eta, and optionally the dimension d.");
      return "ilike_model_dlkj(arma::mat(" + x + "), " + this->scalar_code(arguments[0], context) + ")";
    }
    if (f == "mvnorm")
    {
      if (arguments.size() != 2)
        this->error(call, "mvnorm takes a mean vector and a covariance matrix.");
      return "dmvnorm(arma::colvec(arma::vectorise(" + x + ")), arma::colvec(arma::vectorise(" + this->matrix_code(arguments[0], context) +
             ")), " + this->matrix_code(arguments[1], context) + ")";
    }
    if ((f == "abc_uniform") || (f == "abc_gaussian"))
    {
      // y ~ abc_uniform(s, epsilon, scale = …, p = 2) / abc_gaussian(s, epsilon, scale = …)
      // (dev/plan.md 3.21 (2)): the ABC kernel K_epsilon(y | s), y observed
      Arguments split = this->split_arguments(call);
      std::set<std::string> allowed = {"summaries", "epsilon", "scale"};
      if (f == "abc_uniform")
        allowed.insert("p");
      this->allow_only(call, split, allowed, 2);
      std::map<std::string, RecipeValue> given;
      std::vector<std::string> order = {"summaries", "epsilon"};
      for (size_t k = 0; k < split.positional.size(); ++k)
        given[order[k]] = split.positional[k];
      for (const auto &kv : split.named)
        given[kv.first] = kv.second;
      if (!given.count("summaries") || !given.count("epsilon"))
        this->error(call, f + "(s, epsilon, …): the simulated data (or summaries) s and the tolerance epsilon are needed.");
      std::string scale = given.count("scale") ? "arma::colvec(arma::vectorise(" + this->matrix_code(given["scale"], context) + "))"
                                               : std::string("arma::colvec(1, arma::fill::ones)");
      std::string head = "(arma::colvec(arma::vectorise(" + x + ")), arma::colvec(arma::vectorise(" +
                         this->matrix_code(given["summaries"], context) + ")), " + this->scalar_code(given["epsilon"], context) + ", " + scale;
      if (f == "abc_gaussian")
        return "ilike_model_dabc_gaussian" + head + ")";
      return "ilike_model_dabc_uniform" + head + ", " + (given.count("p") ? this->scalar_code(given["p"], context) : std::string("2.0")) + ")";
    }
    auto found = densities.find(f);
    if (found == densities.end())
      this->error(call, "unknown distribution " + f + " in a model statement (norm, lnorm, gamma, invgamma, unif, exp, pois, tnorm, beta, t, cauchy, halfnorm, halfcauchy, halft, vonmises, chisq, weibull, laplace, logis, bern, binom, negbinom, geom, categorical, mvnorm, mvt, multinom, dirichlet, wishart, invwishart, lkj).");
    if (arguments.size() != found->second.first)
      this->error(call, f + " takes " + std::to_string(found->second.first) + " argument(s).");
    std::string code = found->second.second + "(" + x;
    for (size_t k = 0; k < arguments.size(); ++k)
      code += ", " + this->scalar_code(arguments[k], context);
    return code + ")";
  }

  // distributions whose target is a whole vector (mvnorm, multinom,
  // dirichlet, mvt) or matrix (wishart, invwishart, lkj), written without an
  // index (C227)
  static bool whole_target(const std::string &f)
  {
    return (f == "mvnorm") || (f == "multinom") || (f == "dirichlet") || (f == "mvt") ||
           (f == "wishart") || (f == "invwishart") || (f == "lkj") || (f == "abc_uniform") || (f == "abc_gaussian");
  }

  // the dimension of an lkj matrix, for its simulation: lkj(eta, d), or the
  // largest index of the variable in the statements (C227)
  std::string lkj_dimension(const ModelStatement &st)
  {
    const RecipeValue &arguments = st.value["arguments"];
    CodeContext context;
    context.prior = true;
    if (arguments.size() == 2)
      return "(size_t)(" + this->scalar_code(arguments[1], context) + ")";
    std::string rows, cols;
    this->variable_extent(st.variable, rows, cols);
    if ((rows == "1") || rows.empty())
      this->error(st.node, "lkj: the dimension of " + st.variable + " is not known from the statements; give it: lkj(eta, d).");
    return "(size_t)(" + rows + ")";
  }

  // a statement of a likelihood factor: an observed variable's "~", or a
  // "factor name += expression" term
  bool likelihood_term(const ModelStatement &st) const
  {
    return (st.op == "~" && st.observed) || st.op == "+=";
  }

  // one term's log value: the log density of a "~", or the expression of a "+="
  // (C235) a term written in C++: factor name += cpp::f(i, …, reads = [...]),
  // calling double f(const Parameters &parameters, const Data &data, int i, …)
  // with the loop indices (or index expressions) given; the reads declared
  // are the term's dependencies (as a statement's indices are)
  static bool cpp_term(const RecipeValue &value)
  {
    return (kind_of(value) == "call") && (value["namespace"].as_string() == "cpp");
  }
  const RecipeValue &cpp_term_reads(const RecipeValue &call)
  {
    const RecipeValue &arguments = call["arguments"];
    std::vector<std::string> keys = arguments.names();
    for (size_t k = 0; k < arguments.size(); ++k)
      if ((k < keys.size()) && (keys[k] == "reads"))
      {
        if (kind_of(arguments[k]) != "vector")
          this->error(arguments[k], "reads is a list of what the term reads: reads = [mu[group[i]], sigma].");
        return arguments[k]["elements"];
      }
    this->error(call, "a term written in C++ declares what it reads: cpp::" + call["function"].as_string() +
                "(i, reads = [mu[group[i]], sigma]) (each a variable, all of it, or an element).");
    return arguments; // (not reached)
  }
  std::string cpp_term_code(const RecipeValue &call, const CodeContext &context)
  {
    const RecipeValue &arguments = call["arguments"];
    std::vector<std::string> keys = arguments.names();
    // (::f: the user's function is at global scope, and the generated function
    // around this call has arguments, such as term, that would hide it)
    std::string code = "::" + call["function"].as_string() + "(parameters, " + context.data;
    for (size_t k = 0; k < arguments.size(); ++k)
    {
      std::string key = (k < keys.size()) ? keys[k] : "";
      if (key == "reads")
        continue;
      if (!key.empty())
        this->error(arguments[k], "a term written in C++ takes index arguments and reads = [...] only.");
      code += ", (int)std::lround(" + this->scalar_code(arguments[k], context) + ")";
    }
    this->cpp_term_reads(call); // (checked)
    return "(double)" + code + ")";
  }

  std::string term_code(const ModelStatement &st, const CodeContext &context)
  {
    if ((st.op == "+=") && cpp_term(st.value))
      return this->cpp_term_code(st.value, context);
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
      {"bern", "ilike_model_rbern"}, {"binom", "ilike_model_rbinom"},
      {"invgamma", "ilike_model_rinvgamma"}, {"t", "ilike_model_rt"}, {"cauchy", "ilike_model_rcauchy"},
      {"halfnorm", "ilike_model_rhalfnorm"}, {"halfcauchy", "ilike_model_rhalfcauchy"}, {"halft", "ilike_model_rhalft"},
      {"negbinom", "ilike_model_rnbinom"}, {"geom", "ilike_model_rgeom"}, {"vonmises", "ilike_model_rvonmises"},
      {"chisq", "ilike_model_rchisq"}, {"weibull", "ilike_model_rweibull"}, {"laplace", "ilike_model_rlaplace"},
      {"logis", "ilike_model_rlogis"}};
    auto vector_argument = [&](size_t k) {
      return "arma::colvec(arma::vectorise(" + this->matrix_code(arguments[k], context) + "))"; };
    if (f == "categorical")
      return "ilike_model_rcategorical(rng, " + vector_argument(0) + ")";
    if (f == "multinom")
      return "ilike_model_rmultinom(rng, " + this->scalar_code(arguments[0], context) + ", " + vector_argument(1) + ")";
    if (f == "dirichlet")
      return "ilike_model_rdirichlet(rng, " + vector_argument(0) + ")";
    if (f == "mvt")
      return "ilike_model_rmvt(rng, " + this->scalar_code(arguments[0], context) + ", " + vector_argument(1) + ", " +
             this->matrix_code(arguments[2], context) + ")";
    if ((f == "wishart") || (f == "invwishart"))
      return "ilike_model_r" + f + "(rng, " + this->scalar_code(arguments[0], context) + ", " + this->matrix_code(arguments[1], context) + ")";
    if (f == "lkj")
      return "ilike_model_rlkj(rng, " + this->scalar_code(arguments[0], context) + ", " + this->lkj_dimension(st) + ")";
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
    if (!writable)
      base = this->hoisted(base, "const arma::mat &", context);
    if (whole_target(st.value["function"].as_string()))
    {
      if (st.indices.size() != 0)
        this->error(st.node, st.value["function"].as_string() + "'s target is a whole vector or matrix (no index).");
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

  std::string node_gradient_call_checked(const std::string &name) const
  {
    if (this->is_cpp_node(name))
      this->error(this->deterministic_statement(name).node, "the gradient through the C++ node " + name + " is not available.");
    return "";
  }

  std::string node_gradient_call(const std::string &name, const RecipeValue &indices, const std::string &adjoint,
                                 const CodeContext &context)
  {
    this->node_gradient_call_checked(name);
    std::string data = context.prior ? std::string("nullptr") : (context.node ? std::string("data") : std::string("&data"));
    std::string call = "ilike_model_node_gradient_" + name + "(ilike_G, " + adjoint + ", parameters, " + data;
    for (size_t k = 0; k < indices.size(); ++k)
      call += ", (int)(" + this->scalar_code(indices[k], context) + ")";
    return call + ");\n";
  }

  // (known: the name of a C++ variable already holding the expression's
  // value, or ""; some functions' derivatives are written from their value,
  // which then saves computing their argument)
  std::string gradient_code(const RecipeValue &expression, const std::string &adjoint, const CodeContext &context,
                            const std::string &known = "")
  {
    if (this->parameters_read(expression, context).empty())
      return "";
    std::string kind = kind_of(expression);
    if (kind == "name")
    {
      const std::string name = expression["name"].as_string();
      if (this->is_deterministic(name))
        return this->node_gradient_call(name, RecipeValue::list(), adjoint, context);
      return "if (arma::mat *ilike_gv = " + this->hoisted("ilike_G.find(" + quoted(name) + ")", "arma::mat *const", context) +
             ") (*ilike_gv)(0) += " + adjoint + ";\n";
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
      return "if (arma::mat *ilike_gv = " + this->hoisted("ilike_G.find(" + quoted(name) + ")", "arma::mat *const", context) +
             ") ilike_model_add(*ilike_gv, " + element + ", " + adjoint + ", " + where + ");\n";
    }
    if (kind == "sum")
    {
      // (C233) each term gets the sum's adjoint
      CodeContext inner;
      std::string loop = this->sum_loop(expression, context, inner);
      std::string t = "ilike_a" + std::to_string(++this->gradient_temporaries);
      return "{\n  const double " + t + " = " + adjoint + ";\n  " + loop + "\n  {\n" +
             indented_code(indented_code(this->gradient_code(expression["body"], t, inner))) + "  }\n}\n";
    }
    std::string t = "ilike_a" + std::to_string(++this->gradient_temporaries);
    std::string open = "{\n  const double " + t + " = " + adjoint + ";\n";
    if (kind == "unary")
      return open + indented_code(this->gradient_code(expression["operand"], "-" + t, context)) + "}\n";
    if (kind == "binary")
    {
      std::string op = expression["operator"].as_string();
      const RecipeValue &left = expression["left"], &right = expression["right"];
      // (each operand computed once, into a temporary, if the derivatives use it)
      std::string a = "ilike_l" + std::to_string(++this->gradient_temporaries);
      std::string b = "ilike_r" + std::to_string(this->gradient_temporaries);
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
      std::string operands;
      if (inner.find(a) != std::string::npos)
        operands += "  const double " + a + " = " + this->scalar_code(left, context) + ";\n";
      if (inner.find(b) != std::string::npos)
        operands += "  const double " + b + " = " + this->scalar_code(right, context) + ";\n";
      return open + operands + indented_code(inner) + "}\n";
    }
    if (kind == "call")
    {
      std::string f = expression["function"].as_string();
      const RecipeValue &arguments = expression["arguments"];
      // (C234) piecewise constant: no contribution (as an index reading a parameter)
      if ((f == "equals") || (f == "step"))
        return "";
      // derivatives written from the function's value K, when it is known
      static const std::map<std::string, std::string> from_value = {
        {"exp", "(K)"}, {"expm1", "((K) + 1.0)"}, {"inv_logit", "(K)*(1.0 - (K))"}, {"ilogit", "(K)*(1.0 - (K))"},
        {"tanh", "(1.0 - (K)*(K))"}, {"sqrt", "(0.5/(K))"}};
      auto by_value = from_value.find(f);
      if (!known.empty() && (by_value != from_value.end()))
      {
        std::string d = by_value->second;
        for (size_t at = d.find('K'); at != std::string::npos; at = d.find('K', at + known.size()))
          d.replace(at, 1, known);
        return open + indented_code(this->gradient_code(arguments[0], t + "*" + d, context)) + "}\n";
      }
      // (the argument computed once, into a temporary)
      std::string u = "ilike_u" + std::to_string(++this->gradient_temporaries);
      std::string argument = "  const double " + u + " = " + this->scalar_code(arguments[0], context) + ";\n";
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
        return open + argument + indented_code(this->gradient_code(arguments[0], t + "*" + d, context)) + "}\n";
      }
      std::string v = this->scalar_code(arguments[1], context);
      if (f == "pow")
        return open + argument + indented_code(this->gradient_code(arguments[0], t + "*" + v + "*std::pow(" + u + ", " + v + " - 1.0)", context) +
                                    this->gradient_code(arguments[1], t + "*std::log(" + u + ")*std::pow(" + u + ", " + v + ")", context)) + "}\n";
      if (f == "min" || f == "max")
      {
        // (the derivative of the argument chosen)
        std::string first = (f == "min") ? "((double)" + u + " <= (double)" + v + ")" : "((double)" + u + " >= (double)" + v + ")";
        return open + argument + "  if " + first + "\n  {\n" + indented_code(indented_code(this->gradient_code(arguments[0], t, context))) +
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
      {
        // (C232) d/dS = -S^-1/2 + S^-1 r r' S^-1/2 (= -S^-1/2 + r r'/2 in terms of r = S^-1 (x - m))
        code += std::string("  if (arma::mat *ilike_gv = ilike_G.find(") + quoted(arguments[1]["name"].as_string()) + ")) "
                "*ilike_gv += -0.5*arma::inv_sympd(ilike_S) + 0.5*ilike_r*ilike_r.t();\n";
      }
      return code + "}\n";
    }
    if ((f == "wishart") || (f == "invwishart") || (f == "lkj"))
    {
      // (C232) the gradient with respect to the matrix (each element treated
      // as free: the transforms' chain rule takes the symmetry into account)
      for (size_t k = 0; k < arguments.size(); ++k)
      {
        std::set<std::string> found;
        this->read_parameters(arguments[k], found);
        if (!found.empty())
          this->error(call, "the gradient with respect to the arguments of " + f + " is not available.");
      }
      if (st.observed)
        return "";
      std::string target = this->target_code(st, context, false);
      std::string code = "{\n  arma::mat ilike_X = " + target + ";\n  arma::mat ilike_Xi;\n"
                         "  if (!arma::inv(ilike_Xi, ilike_X)) ilike_Xi = arma::mat(arma::size(ilike_X), arma::fill::value(arma::datum::nan));\n"
                         "  const double ilike_p = (double)ilike_X.n_rows;\n";
      if (f == "wishart")
        code += "  arma::mat ilike_d = 0.5*(" + this->scalar_code(arguments[0], context) + " - ilike_p - 1.0)*ilike_Xi - 0.5*arma::inv(" +
                this->matrix_code(arguments[1], context) + ");\n";
      else if (f == "invwishart")
        code += "  arma::mat ilike_d = -0.5*(" + this->scalar_code(arguments[0], context) + " + ilike_p + 1.0)*ilike_Xi + 0.5*ilike_Xi*" +
                this->matrix_code(arguments[1], context) + "*ilike_Xi;\n";
      else
        code += "  arma::mat ilike_d = (" + this->scalar_code(arguments[0], context) + " - 1.0)*ilike_Xi;\n";
      code += std::string("  if (arma::mat *ilike_gv = ilike_G.find(") + quoted(st.variable) + ")) *ilike_gv += ilike_d;\n";
      return code + "}\n";
    }
    // (C227) the vector distributions: the gradient with respect to the
    // target (if continuous and a parameter) and to a vector argument that
    // is a parameter written by its name; other arguments that read
    // parameters have no gradient (the factor then has none)
    if ((f == "mvt") || (f == "dirichlet") || (f == "multinom") || (f == "categorical"))
    {
      auto parameter_name = [&](const RecipeValue &e) -> std::string {
        if ((kind_of(e) == "name") && !this->data_names.count(e["name"].as_string()) &&
            !this->constant_names.count(e["name"].as_string()) && !this->is_deterministic(e["name"].as_string()))
          return e["name"].as_string();
        return std::string();
      };
      auto refuse_reading = [&](const RecipeValue &e, const std::string &what) {
        std::set<std::string> found;
        this->read_parameters(e, found);
        if (!found.empty())
          this->error(call, "the gradient with respect to " + what + " of " + f + " is available only for a parameter written by its name.");
      };
      std::string target = this->target_code(st, context, false);
      auto add_to = [&](const std::string &name, const std::string &value) {
        return "  if (arma::mat *ilike_gv = ilike_G.find(" + quoted(name) + ")) *ilike_gv += arma::reshape(" + value +
               ", ilike_gv->n_rows, ilike_gv->n_cols);\n";
      };
      std::string code = "{\n";
      if (f == "mvt")
      {
        code += "  arma::colvec ilike_x = arma::vectorise(" + target + ");\n"
                "  const double ilike_nu = " + this->scalar_code(arguments[0], context) + ";\n"
                "  arma::colvec ilike_m = arma::vectorise(" + this->matrix_code(arguments[1], context) + ");\n"
                "  arma::mat ilike_S = " + this->matrix_code(arguments[2], context) + ";\n"
                "  arma::colvec ilike_r = arma::solve(ilike_S, ilike_x - ilike_m, arma::solve_opts::likely_sympd);\n"
                "  const double ilike_q = arma::dot(ilike_x - ilike_m, ilike_r);\n"
                "  arma::colvec ilike_w = ((ilike_nu + ilike_x.n_elem)/(ilike_nu + ilike_q))*ilike_r;\n";
        if (!st.observed)
          code += "  if (arma::mat *ilike_gv = ilike_G.find(" + quoted(st.variable) + ")) *ilike_gv -= arma::reshape(ilike_w, ilike_gv->n_rows, ilike_gv->n_cols);\n";
        std::string mean = parameter_name(arguments[1]);
        if (!mean.empty())
          code += add_to(mean, "ilike_w");
        else
          refuse_reading(arguments[1], "the location");
        refuse_reading(arguments[0], "the degrees of freedom");
        const std::string scale = parameter_name(arguments[2]);
        if (!scale.empty())
          // (C232) d/dS = -S^-1/2 + (nu + p)/2 S^-1 r r' S^-1/(nu + q), r = x - m
          code += std::string("  if (arma::mat *ilike_gv = ilike_G.find(") + quoted(scale) + ")) *ilike_gv += -0.5*arma::inv_sympd(ilike_S) + "
                  "0.5*((ilike_nu + ilike_x.n_elem)/(ilike_nu + ilike_q))*ilike_r*ilike_r.t();\n";
        else
          refuse_reading(arguments[2], "the scale matrix");
      }
      else if (f == "dirichlet")
      {
        code += "  arma::colvec ilike_x = arma::vectorise(" + target + ");\n"
                "  arma::colvec ilike_a = arma::vectorise(" + this->matrix_code(arguments[0], context) + ");\n";
        if (!st.observed)
          code += add_to(st.variable, "arma::colvec((ilike_a - 1.0)/ilike_x)");
        std::string alpha = parameter_name(arguments[0]);
        if (!alpha.empty())
          code += "  {\n    arma::colvec ilike_d(ilike_a.n_elem);\n    const double ilike_s = R::digamma(arma::accu(ilike_a));\n"
                  "    for (arma::uword k = 0; k < ilike_a.n_elem; ++k) ilike_d[k] = ilike_s - R::digamma(ilike_a[k]) + std::log(ilike_x[k]);\n  " +
                  add_to(alpha, "ilike_d") + "  }\n";
        else
          refuse_reading(arguments[0], "the concentrations");
      }
      else if (f == "multinom")
      {
        // (the target is discrete)
        code += "  arma::colvec ilike_x = arma::vectorise(" + target + ");\n"
                "  arma::colvec ilike_p = arma::vectorise(" + this->matrix_code(arguments[1], context) + ");\n";
        std::string probabilities = parameter_name(arguments[1]);
        if (!probabilities.empty())
          code += add_to(probabilities, "arma::colvec(ilike_x/ilike_p)");
        else
          refuse_reading(arguments[1], "the probabilities");
      }
      else
      {
        // categorical (the target is discrete): d log p[x] / d p[x] = 1/p[x]
        code += "  const double ilike_x = " + target + ";\n"
                "  arma::colvec ilike_p = arma::vectorise(" + this->matrix_code(arguments[0], context) + ");\n";
        std::string probabilities = parameter_name(arguments[0]);
        if (!probabilities.empty())
          code += "  {\n    arma::colvec ilike_d(ilike_p.n_elem, arma::fill::zeros);\n"
                  "    arma::uword ilike_k = (arma::uword)(ilike_x - 1.0);\n"
                  "    if (ilike_k < ilike_p.n_elem) ilike_d[ilike_k] = 1.0/ilike_p[ilike_k];\n  " + add_to(probabilities, "ilike_d") + "  }\n";
        else
          refuse_reading(arguments[0], "the probabilities");
      }
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
    else if (f == "invgamma")
    {
      // log f = a log b - lgamma(a) - (a + 1) log x - b/x
      dx = "(-(" + v[0] + " + 1.0)/" + x + " + " + v[1] + "/(" + x + "*" + x + "))";
      d[0] = "(std::log(" + v[1] + ") - R::digamma(" + v[0] + ") - std::log(" + x + "))";
      d[1] = "(" + v[0] + "/" + v[1] + " - 1.0/" + x + ")";
    }
    else if ((f == "t") || (f == "halft"))
    {
      // t(nu, m, s) (halft(nu, s): m = 0, x >= 0): z = (x - m)/s,
      // log f = lgamma((nu+1)/2) - lgamma(nu/2) - log(nu pi)/2 - log s - (nu+1)/2 log(1 + z^2/nu)
      std::string m = (f == "t") ? v[1] : "0.0", sd = (f == "t") ? v[2] : v[1];
      prelude = "  const double ilike_z = (" + x + " - " + m + ")/" + sd + ";\n";
      std::string q = "(" + v[0] + " + ilike_z*ilike_z)";
      dx = "(-(" + v[0] + " + 1.0)*ilike_z/(" + sd + "*" + q + "))";
      d[0] = "(0.5*(R::digamma((" + v[0] + " + 1.0)/2.0) - R::digamma(" + v[0] + "/2.0) - 1.0/" + v[0] +
             " - std::log1p(ilike_z*ilike_z/" + v[0] + ") + (" + v[0] + " + 1.0)*ilike_z*ilike_z/(" + v[0] + "*" + q + ")))";
      std::string dsd = "(-1.0/" + sd + " + (" + v[0] + " + 1.0)*ilike_z*ilike_z/(" + sd + "*" + q + "))";
      if (f == "t")
      {
        d[1] = "(-" + dx + ")";
        d[2] = dsd;
      }
      else
        d[1] = dsd;
    }
    else if ((f == "cauchy") || (f == "halfcauchy"))
    {
      // cauchy(m, s) (halfcauchy(s): m = 0, x >= 0): log f = -log(pi s) - log(1 + z^2)
      std::string m = (f == "cauchy") ? v[0] : "0.0", sd = (f == "cauchy") ? v[1] : v[0];
      prelude = "  const double ilike_z = (" + x + " - " + m + ")/" + sd + ";\n";
      dx = "(-2.0*ilike_z/(" + sd + "*(1.0 + ilike_z*ilike_z)))";
      std::string dsd = "(-1.0/" + sd + " + 2.0*ilike_z*ilike_z/(" + sd + "*(1.0 + ilike_z*ilike_z)))";
      if (f == "cauchy")
      {
        d[0] = "(-" + dx + ")";
        d[1] = dsd;
      }
      else
        d[0] = dsd;
    }
    else if (f == "halfnorm")
    {
      dx = "(-" + x + "/(" + v[0] + "*" + v[0] + "))";
      d[0] = "(-1.0/" + v[0] + " + " + x + "*" + x + "/(" + v[0] + "*" + v[0] + "*" + v[0] + "))";
    }
    else if (f == "negbinom")
    {
      // (the target is discrete) log f = lgamma(x + n) - lgamma(n) - lgamma(x + 1) + n log p + x log(1 - p)
      d[0] = "(R::digamma(" + x + " + " + v[0] + ") - R::digamma(" + v[0] + ") + std::log(" + v[1] + "))";
      d[1] = "(" + v[0] + "/" + v[1] + " - " + x + "/(1.0 - " + v[1] + "))";
    }
    else if (f == "geom")
      d[0] = "(1.0/" + v[0] + " - " + x + "/(1.0 - " + v[0] + "))"; // (the target is discrete)
    else if (f == "vonmises")
    {
      // log f = kappa cos(x - mu) - log(2 pi I0(kappa))
      dx = "(-" + v[1] + "*std::sin(" + x + " - " + v[0] + "))";
      d[0] = "(" + v[1] + "*std::sin(" + x + " - " + v[0] + "))";
      d[1] = "(std::cos(" + x + " - " + v[0] + ") - R::bessel_i(" + v[1] + ", 1.0, 2.0)/R::bessel_i(" + v[1] + ", 0.0, 2.0))";
    }
    else if (f == "chisq")
    {
      // (C229) log f = (k/2 - 1) log x - x/2 - (k/2) log 2 - lgamma(k/2)
      dx = "((" + v[0] + "/2.0 - 1.0)/" + x + " - 0.5)";
      d[0] = "(0.5*(std::log(" + x + ") - std::log(2.0) - R::digamma(" + v[0] + "/2.0)))";
    }
    else if (f == "weibull")
    {
      // log f = log k - log l + (k - 1) log(x/l) - (x/l)^k
      prelude = "  const double ilike_z = " + x + "/" + v[1] + ";\n  const double ilike_zk = std::pow(ilike_z, " + v[0] + ");\n";
      dx = "((" + v[0] + " - 1.0)/" + x + " - " + v[0] + "*ilike_zk/" + x + ")";
      d[0] = "(1.0/" + v[0] + " + std::log(ilike_z) - ilike_zk*std::log(ilike_z))";
      d[1] = "((-" + v[0] + " + " + v[0] + "*ilike_zk)/" + v[1] + ")";
    }
    else if (f == "laplace")
    {
      // log f = -log(2 b) - |x - m|/b
      std::string sign = "((" + x + " > " + v[0] + ") ? 1.0 : ((" + x + " < " + v[0] + ") ? -1.0 : 0.0))";
      dx = "(-" + sign + "/" + v[1] + ")";
      d[0] = "(" + sign + "/" + v[1] + ")";
      d[1] = "(-1.0/" + v[1] + " + std::abs(" + x + " - " + v[0] + ")/(" + v[1] + "*" + v[1] + "))";
    }
    else if (f == "logis")
    {
      // d log f / d z = -tanh(z/2), z = (x - m)/s
      prelude = "  const double ilike_z = (" + x + " - " + v[0] + ")/" + v[1] + ";\n  const double ilike_t = std::tanh(0.5*ilike_z);\n";
      dx = "(-ilike_t/" + v[1] + ")";
      d[0] = "(ilike_t/" + v[1] + ")";
      d[1] = "(-1.0/" + v[1] + " + ilike_z*ilike_t/" + v[1] + ")";
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
        code += indented_code(this->gradient_code(arguments[k], d[k], context, (k < 4) ? v[k] : std::string()));
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
      // the log density and the gradients in one pass over the terms (C222):
      // every term's value, and the gradients of those reading a variable asked for
      std::string value_body = "ilike_model_gradients ilike_G(variables, gradients, parameters);\ndouble ilike_value = 0.0;\ntry\n{\n";
      for (const ModelStatement* st : statements)
      {
        std::set<std::string> reads = this->statements_read(std::vector<const ModelStatement*>(1, st));
        CodeContext context;
        context.prior = prior;
        for (const ModelLoop &loop : st->loops) context.indices.insert(loop.index);
        // (the lookups by name of the value and of the gradient, made once
        // before the statement's loops; separately for the two functions, whose
        // code is in different places)
        Hoist value_hoist, gradient_hoist;
        CodeContext value_context = context, gradient_context = context;
        value_context.hoist = &value_hoist;
        gradient_context.hoist = &gradient_hoist;
        std::string value = "ilike_value += " + this->term_code(*st, value_context) + ";\n";
        CodeContext fresh;
        fresh.prior = prior;
        if (reads.empty())
        {
          value_body += indented_code("{\n" + value_hoist.declarations + this->loops_code(*st, fresh, value, "") + "}\n");
          continue;
        }
        std::string guard;
        for (const std::string &r : reads)
          guard += (guard.empty() ? "" : " || ") + std::string("ilike_G.find(") + quoted(r) + ")";
        std::string term = this->term_gradient_code(*st, gradient_context);
        body += "  if (" + guard + ")\n  {\n" +
                indented_code(indented_code(gradient_hoist.declarations + this->loops_code(*st, fresh, term, ""))) + "  }\n";
        std::string value_term = this->term_gradient_code(*st, value_context);
        CodeContext fresh2;
        fresh2.prior = prior;
        CodeContext fresh3;
        fresh3.prior = prior;
        value_body += "  {\n" + indented_code(indented_code(value_hoist.declarations)) +
                      "  if (" + guard + ")\n  {\n" + indented_code(indented_code(this->loops_code(*st, fresh2, value + value_term, ""))) +
                      "  }\n  else\n  {\n" + indented_code(indented_code(this->loops_code(*st, fresh3, value, ""))) + "  }\n  }\n";
      }
      body += "}\ncatch (const ilike_model_outside &)\n{\n  for (arma::mat &g : gradients) g.fill(arma::datum::nan); // (density 0 at this point)\n}\n";
      value_body += "}\ncatch (const ilike_model_outside &)\n{\n  for (arma::mat &g : gradients) g.fill(arma::datum::nan);\n"
                    "  return -std::numeric_limits<double>::infinity(); // (density 0 at this point)\n}\nreturn ilike_value;\n";
      std::string many = "(const std::vector<std::string> &variables, std::vector<arma::mat> &gradients, const Parameters &parameters" +
                         std::string(prior ? ")" : ", const Data &data)");
      std::string code = "void " + name + "_many" + many + "\n{\n" + indented_code(body) + "}\n\n";
      code += "double " + name + "_many_value" + many + "\n{\n" + indented_code(value_body) + "}\n\n";
      std::string signature = prior ? "(const std::string &variable, const Parameters &parameters)"
                                    : "(const std::string &variable, const Parameters &parameters, const Data &data)";
      code += "arma::mat " + name + signature + "\n{\n  std::vector<std::string> variables(1, variable);\n  std::vector<arma::mat> gradients;\n  " +
              name + "_many(variables, gradients, parameters" + (prior ? "" : ", data") + ");\n  return gradients[0];\n}\n\n";
      if (prior)
        code += "arma::mat " + name + "_plate(const std::string &variable, const Parameters &parameters, const Data &data)\n{\n"
                "  return " + name + "(variable, parameters);\n}\n\n"
                "void " + name + "_many_plate(const std::vector<std::string> &variables, std::vector<arma::mat> &gradients, "
                "const Parameters &parameters, const Data &data)\n{\n  " + name + "_many(variables, gradients, parameters);\n}\n\n"
                "double " + name + "_many_value_plate(const std::vector<std::string> &variables, std::vector<arma::mat> &gradients, "
                "const Parameters &parameters, const Data &data)\n{\n  return " + name + "_many_value(variables, gradients, parameters);\n}\n\n";
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
    if (kind == "sum")
    {
      // (C233) the reads of every term
      CodeContext inner;
      std::string loop = this->sum_loop(expression, context, inner);
      std::string body = this->reads_code(expression["body"], inner, whole);
      return body.empty() ? std::string() : loop + "\n{\n" + indented_code(body) + "}\n";
    }
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
    if ((st.op == "+=") && cpp_term(st.value))
    {
      // (C235) the declared reads: a variable (all of it) or an element
      std::string code;
      const RecipeValue &reads = this->cpp_term_reads(st.value);
      for (size_t k = 0; k < reads.size(); ++k)
      {
        std::string kind = kind_of(reads[k]);
        if ((kind != "name") && (kind != "index"))
          this->error(reads[k], "each of reads = [...] is a variable or an element (x or x[i]).");
        code += this->reads_code(reads[k], context, true);
      }
      return code;
    }
    if (st.op == "+=")
      return this->reads_code(st.value, context, false);
    std::string code;
    bool mvnorm = whole_target(st.value["function"].as_string());
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

      Hoist hoist;
      CodeContext term_context = context;
      term_context.hoist = &hoist;
      std::string term_value = this->term_code(st, term_context);
      term += "    case " + std::to_string(s) + ":\n    {\n" + indented_code(indented_code(indented_code(hoist.declarations)));
      for (size_t k = 0; k < st.loops.size(); ++k)
        term += "      int " + st.loops[k].index + " = (int)terms(term, " + std::to_string(k + 1) + ");\n";
      term += "      return " + term_value + ";\n    }\n";
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
          std::string assign = whole_target(st.value["function"].as_string())
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
// the partial derivatives of log dtnorm(x, m, s, l, u) with respect to x, m, s, l, u
inline arma::colvec ilike_model_dtnorm_partials(double x, double m, double s, double l, double u)
{
  // phi(alpha)/Z and phi(beta)/Z as exp(log phi - log Z), with log Z from
  // log_normal_mass: accurate when both limits are far out on one side (C218;
  // Z itself would underflow or cancel there)
  double alpha = (l - m)/s, beta = (u - m)/s;
  double log_Z = log_normal_mass(alpha, beta);
  double ra = std::isfinite(alpha) ? std::exp(R::dnorm(alpha, 0.0, 1.0, 1) - log_Z) : 0.0;
  double rb = std::isfinite(beta) ? std::exp(R::dnorm(beta, 0.0, 1.0, 1) - log_Z) : 0.0;
  double a_ra = std::isfinite(alpha) ? alpha*ra : 0.0, b_rb = std::isfinite(beta) ? beta*rb : 0.0;
  arma::colvec d(5);
  d[0] = -(x - m)/(s*s);
  d[1] = (x - m)/(s*s) + (rb - ra)/s;
  d[2] = -1.0/s + (x - m)*(x - m)/(s*s*s) + (b_rb - a_ra)/s;
  d[3] = ra/s;
  d[4] = -rb/s;
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
// (C214) invgamma(a, b): density b^a / Gamma(a) x^(-a-1) exp(-b/x), x > 0
// (1/x ~ gamma(a, rate b)); t(nu, m, s): Student t, location m, scale s;
// cauchy(m, s); halfnorm(s), halfcauchy(s), halft(nu, s): twice the
// density of the zero-centred distribution on x >= 0
inline double ilike_model_dinvgamma(double x, double a, double b)
{
  if (x <= 0.0 || a <= 0.0 || b <= 0.0) return -std::numeric_limits<double>::infinity();
  return a*std::log(b) - std::lgamma(a) - (a + 1.0)*std::log(x) - b/x;
}
inline double ilike_model_dt(double x, double nu, double m, double s)
{
  if (nu <= 0.0 || s <= 0.0) return -std::numeric_limits<double>::infinity();
  double z = (x - m)/s;
  return std::lgamma((nu + 1.0)/2.0) - std::lgamma(nu/2.0) - 0.5*std::log(nu*3.14159265358979323846) - std::log(s) -
         0.5*(nu + 1.0)*std::log1p(z*z/nu);
}
inline double ilike_model_dcauchy(double x, double m, double s)
{
  if (s <= 0.0) return -std::numeric_limits<double>::infinity();
  double z = (x - m)/s;
  return -std::log(3.14159265358979323846*s) - std::log1p(z*z);
}
inline double ilike_model_dhalfnorm(double x, double s)
{
  if (x < 0.0 || s <= 0.0) return -std::numeric_limits<double>::infinity();
  return std::log(2.0) + dnorm(x, 0.0, s);
}
inline double ilike_model_dhalfcauchy(double x, double s)
{
  if (x < 0.0) return -std::numeric_limits<double>::infinity();
  return std::log(2.0) + ilike_model_dcauchy(x, 0.0, s);
}
inline double ilike_model_dhalft(double x, double nu, double s)
{
  if (x < 0.0) return -std::numeric_limits<double>::infinity();
  return std::log(2.0) + ilike_model_dt(x, nu, 0.0, s);
}
inline double ilike_model_rinvgamma(RandomNumberGenerator &rng, double a, double b)
{
  return 1.0/rgamma(rng, a, b);
}
inline double ilike_model_rt(RandomNumberGenerator &rng, double nu, double m, double s)
{
  // a normal over the square root of a chi-squared on nu degrees of freedom, divided by nu
  double z = rnorm(rng, 0.0, 1.0), w = rgamma(rng, nu/2.0, 0.5);
  return m + s*z/std::sqrt(w/nu);
}
inline double ilike_model_rcauchy(RandomNumberGenerator &rng, double m, double s)
{
  return m + s*std::tan(3.14159265358979323846*(runif(rng) - 0.5));
}
inline double ilike_model_rhalfnorm(RandomNumberGenerator &rng, double s)
{
  return std::abs(rnorm(rng, 0.0, s));
}
inline double ilike_model_rhalfcauchy(RandomNumberGenerator &rng, double s)
{
  return std::abs(ilike_model_rcauchy(rng, 0.0, s));
}
inline double ilike_model_rhalft(RandomNumberGenerator &rng, double nu, double s)
{
  return std::abs(ilike_model_rt(rng, nu, 0.0, s));
}
// (C227) negbinom(size, prob) (R's), geom(prob) (failures before the first
// success, R's), categorical(p) (x in 1, ..., K, p the probabilities),
// vonmises(mu, kappa); the vector and matrix distributions multinom(n, p),
// dirichlet(alpha), mvt(df, location, scale), wishart(df, S),
// invwishart(df, S) (R's parameterisations, as MCMCpack's riwish) and
// lkj(eta) (correlation matrices, Lewandowski, Kurowicka & Joe 2009)
// (C229) chisq(df), weibull(shape, scale), logis(location, scale) (R's);
// laplace(location, scale)
inline double ilike_model_dchisq(double x, double k)
{
  if (!(x > 0.0) || !(k > 0.0)) return -std::numeric_limits<double>::infinity();
  return (0.5*k - 1.0)*std::log(x) - 0.5*x - 0.5*k*std::log(2.0) - std::lgamma(0.5*k);
}
inline double ilike_model_rchisq(RandomNumberGenerator &rng, double k)
{
  return rgamma(rng, 0.5*k, 0.5);
}
inline double ilike_model_dweibull(double x, double k, double l)
{
  if (!(x > 0.0) || !(k > 0.0) || !(l > 0.0)) return -std::numeric_limits<double>::infinity();
  double z = x/l;
  return std::log(k/l) + (k - 1.0)*std::log(z) - std::pow(z, k);
}
inline double ilike_model_rweibull(RandomNumberGenerator &rng, double k, double l)
{
  return l*std::pow(-std::log(runif(rng)), 1.0/k);
}
inline double ilike_model_dlaplace(double x, double m, double b)
{
  if (!(b > 0.0)) return -std::numeric_limits<double>::infinity();
  return -std::log(2.0*b) - std::abs(x - m)/b;
}
inline double ilike_model_rlaplace(RandomNumberGenerator &rng, double m, double b)
{
  // the difference of two exponentials
  return m + b*(std::log(runif(rng)) - std::log(runif(rng)));
}
inline double ilike_model_dlogis(double x, double m, double s)
{
  if (!(s > 0.0)) return -std::numeric_limits<double>::infinity();
  double z = std::abs((x - m)/s);
  return -std::log(s) - z - 2.0*std::log1p(std::exp(-z));
}
inline double ilike_model_rlogis(RandomNumberGenerator &rng, double m, double s)
{
  double u = runif(rng);
  return m + s*(std::log(u) - std::log1p(-u));
}
inline double ilike_model_dnbinom(double x, double n, double p)
{
  if (x < 0.0 || x != std::floor(x) || !(n > 0.0) || !(p > 0.0) || p > 1.0) return -std::numeric_limits<double>::infinity();
  return std::lgamma(x + n) - std::lgamma(n) - std::lgamma(x + 1.0) + n*std::log(p) + ((x > 0.0) ? x*std::log1p(-p) : 0.0);
}
inline double ilike_model_rnbinom(RandomNumberGenerator &rng, double n, double p)
{
  // a Poisson with a gamma-distributed mean
  if (p >= 1.0) return 0.0;
  return (double)rpois(rng, rgamma(rng, n, p/(1.0 - p)));
}
inline double ilike_model_dgeom(double x, double p)
{
  if (x < 0.0 || x != std::floor(x) || !(p > 0.0) || p > 1.0) return -std::numeric_limits<double>::infinity();
  return std::log(p) + ((x > 0.0) ? x*std::log1p(-p) : 0.0);
}
inline double ilike_model_rgeom(RandomNumberGenerator &rng, double p)
{
  if (p >= 1.0) return 0.0;
  return std::floor(std::log(runif(rng))/std::log1p(-p));
}
inline double ilike_model_dcategorical(double x, const arma::colvec &p)
{
  if (x < 1.0 || x > double(p.n_elem) || x != std::floor(x)) return -std::numeric_limits<double>::infinity();
  return std::log(p[(arma::uword)(x - 1.0)]);
}
inline double ilike_model_rcategorical(RandomNumberGenerator &rng, const arma::colvec &p)
{
  double u = runif(rng)*arma::accu(p), cumulative = 0.0;
  for (arma::uword k = 0; k < p.n_elem; ++k)
  {
    cumulative = cumulative + p[k];
    if (u < cumulative) return double(k + 1);
  }
  return double(p.n_elem);
}
inline double ilike_model_dvonmises(double x, double mu, double kappa)
{
  if (!(kappa >= 0.0)) return -std::numeric_limits<double>::infinity();
  // log I0(kappa) from the exponentially scaled Bessel function
  return kappa*(std::cos(x - mu) - 1.0) - std::log(2.0*3.14159265358979323846*R::bessel_i(kappa, 0.0, 2.0));
}
inline double ilike_model_rvonmises(RandomNumberGenerator &rng, double mu, double kappa)
{
  // Best & Fisher (1979)
  const double pi = 3.14159265358979323846;
  if (kappa < 1e-8) return mu + pi*(2.0*runif(rng) - 1.0);
  double tau = 1.0 + std::sqrt(1.0 + 4.0*kappa*kappa);
  double rho = (tau - std::sqrt(2.0*tau))/(2.0*kappa);
  double r = (1.0 + rho*rho)/(2.0*rho);
  double f;
  while (true)
  {
    double u1 = runif(rng), u2 = runif(rng);
    double z = std::cos(pi*u1);
    f = (1.0 + r*z)/(r + z);
    double c = kappa*(r - f);
    if ((c*(2.0 - c) - u2 > 0.0) || (std::log(c/u2) + 1.0 - c >= 0.0)) break;
  }
  double u3 = runif(rng);
  return mu + ((u3 > 0.5) ? 1.0 : -1.0)*std::acos(std::max(-1.0, std::min(1.0, f)));
}
inline double ilike_model_dmultinom(const arma::colvec &x, double n, const arma::colvec &p)
{
  if (x.n_elem != p.n_elem) Rcpp::stop("multinom: the counts and the probabilities differ in length.");
  double total = 0.0, value = std::lgamma(n + 1.0);
  for (arma::uword k = 0; k < x.n_elem; ++k)
  {
    if (x[k] < 0.0 || x[k] != std::floor(x[k])) return -std::numeric_limits<double>::infinity();
    total = total + x[k];
    value = value - std::lgamma(x[k] + 1.0) + ((x[k] > 0.0) ? x[k]*std::log(p[k]) : 0.0);
  }
  return (total == n) ? value : -std::numeric_limits<double>::infinity();
}
inline arma::mat ilike_model_rmultinom(RandomNumberGenerator &rng, double n, const arma::colvec &p)
{
  arma::colvec x(p.n_elem, arma::fill::zeros);
  double total = arma::accu(p);
  for (int i = 0; i < (int)n; ++i)
    x[(arma::uword)(ilike_model_rcategorical(rng, p/total) - 1.0)] += 1.0;
  return x;
}
inline double ilike_model_ddirichlet(const arma::colvec &x, const arma::colvec &alpha)
{
  if (x.n_elem != alpha.n_elem) Rcpp::stop("dirichlet: the vector and the concentrations differ in length.");
  if (arma::any(x <= 0.0) || std::abs(arma::accu(x) - 1.0) > 1e-8) return -std::numeric_limits<double>::infinity();
  double value = std::lgamma(arma::accu(alpha));
  for (arma::uword k = 0; k < x.n_elem; ++k)
    value = value - std::lgamma(alpha[k]) + (alpha[k] - 1.0)*std::log(x[k]);
  return value;
}
inline arma::mat ilike_model_rdirichlet(RandomNumberGenerator &rng, const arma::colvec &alpha)
{
  arma::colvec x(alpha.n_elem);
  for (arma::uword k = 0; k < alpha.n_elem; ++k) x[k] = rgamma(rng, alpha[k], 1.0);
  return x/arma::accu(x);
}
inline double ilike_model_dmvt(const arma::colvec &x, double nu, const arma::colvec &m, const arma::mat &S)
{
  double p = double(x.n_elem), log_det, sign;
  arma::log_det(log_det, sign, S);
  if (!(nu > 0.0) || !(sign > 0.0)) return -std::numeric_limits<double>::infinity();
  arma::colvec r = arma::solve(S, x - m, arma::solve_opts::likely_sympd);
  double q = arma::dot(x - m, r);
  return std::lgamma(0.5*(nu + p)) - std::lgamma(0.5*nu) - 0.5*p*std::log(nu*3.14159265358979323846) - 0.5*log_det -
         0.5*(nu + p)*std::log1p(q/nu);
}
inline arma::mat ilike_model_rmvt(RandomNumberGenerator &rng, double nu, const arma::colvec &m, const arma::mat &S)
{
  arma::colvec z(m.n_elem);
  for (arma::uword k = 0; k < m.n_elem; ++k) z[k] = rnorm(rng, 0.0, 1.0);
  double w = rgamma(rng, 0.5*nu, 0.5);
  return m + arma::chol(S, "lower")*z/std::sqrt(w/nu);
}
// log of the multivariate gamma function Gamma_p(a)
inline double ilike_model_lmvgamma(double a, double p)
{
  double value = 0.25*p*(p - 1.0)*std::log(3.14159265358979323846);
  for (int j = 1; j <= (int)p; ++j) value = value + std::lgamma(a + 0.5*(1.0 - j));
  return value;
}
inline bool ilike_model_sympd(const arma::mat &X)
{
  return (X.n_rows == X.n_cols) && X.is_symmetric(1e-10*(1.0 + arma::abs(X).max())) && arma::mat(X).is_sympd();
}
inline double ilike_model_dwishart(const arma::mat &X, double nu, const arma::mat &S)
{
  double p = double(X.n_rows), log_det_X, log_det_S, sign;
  if (!ilike_model_sympd(X) || !(nu > p - 1.0)) return -std::numeric_limits<double>::infinity();
  arma::log_det(log_det_X, sign, X);
  arma::log_det(log_det_S, sign, S);
  return 0.5*(nu - p - 1.0)*log_det_X - 0.5*arma::trace(arma::solve(S, X)) - 0.5*nu*p*std::log(2.0) - 0.5*nu*log_det_S -
         ilike_model_lmvgamma(0.5*nu, p);
}
inline arma::mat ilike_model_rwishart(RandomNumberGenerator &rng, double nu, const arma::mat &S)
{
  // Bartlett: A lower triangular, A_ii^2 ~ chisq(nu - i + 1), A_ij ~ N(0, 1)
  arma::uword p = S.n_rows;
  arma::mat A(p, p, arma::fill::zeros);
  for (arma::uword i = 0; i < p; ++i)
  {
    A(i, i) = std::sqrt(rgamma(rng, 0.5*(nu - double(i)), 0.5));
    for (arma::uword j = 0; j < i; ++j) A(i, j) = rnorm(rng, 0.0, 1.0);
  }
  arma::mat LA = arma::chol(S, "lower")*A;
  return LA*LA.t();
}
inline double ilike_model_dinvwishart(const arma::mat &X, double nu, const arma::mat &Psi)
{
  double p = double(X.n_rows), log_det_X, log_det_Psi, sign;
  if (!ilike_model_sympd(X) || !(nu > p - 1.0)) return -std::numeric_limits<double>::infinity();
  arma::log_det(log_det_X, sign, X);
  arma::log_det(log_det_Psi, sign, Psi);
  return 0.5*nu*log_det_Psi - 0.5*nu*p*std::log(2.0) - ilike_model_lmvgamma(0.5*nu, p) - 0.5*(nu + p + 1.0)*log_det_X -
         0.5*arma::trace(arma::solve(X, Psi));
}
inline arma::mat ilike_model_rinvwishart(RandomNumberGenerator &rng, double nu, const arma::mat &Psi)
{
  arma::mat W = ilike_model_rwishart(rng, nu, arma::inv_sympd(Psi));
  return arma::inv_sympd(W);
}
inline double ilike_model_dlkj(const arma::mat &R, double eta)
{
  // density det(R)^(eta - 1)/c_d(eta), Lewandowski, Kurowicka & Joe (2009):
  // c_d = 2^(sum_k (2 eta - 2 + d - k)(d - k)) prod_k B(b_k, b_k)^(d - k),
  // b_k = eta + (d - k - 1)/2, k = 1, ..., d - 1
  double d = double(R.n_rows), log_det, sign;
  if (!(eta > 0.0) || !ilike_model_sympd(R) || arma::any(arma::abs(R.diag() - 1.0) > 1e-10))
    return -std::numeric_limits<double>::infinity();
  arma::log_det(log_det, sign, R);
  double log_c = 0.0;
  for (int k = 1; k < (int)d; ++k)
  {
    double b = eta + 0.5*(d - k - 1.0);
    log_c = log_c + (2.0*eta - 2.0 + d - k)*(d - k)*std::log(2.0) + (d - k)*(2.0*std::lgamma(b) - std::lgamma(2.0*b));
  }
  return (eta - 1.0)*log_det - log_c;
}
inline arma::mat ilike_model_rlkj(RandomNumberGenerator &rng, double eta, size_t d)
{
  // the C-vine method of Lewandowski, Kurowicka & Joe (2009): partial
  // correlations P[k][i] ~ 2 Beta(b_k, b_k) - 1, b_k = eta + (d - 1 - k)/2
  // (k = 1, ..., d - 1), turned into correlations
  arma::mat P(d, d, arma::fill::zeros), S(d, d, arma::fill::eye);
  for (size_t k = 0; k + 1 < d; ++k)
  {
    double b = eta + 0.5*(double(d) - 2.0 - double(k));
    for (size_t i = k + 1; i < d; ++i)
    {
      double g1 = rgamma(rng, b, 1.0), g2 = rgamma(rng, b, 1.0);
      P(k, i) = 2.0*g1/(g1 + g2) - 1.0;
      double p = P(k, i);
      for (int l = int(k) - 1; l >= 0; --l)
        p = p*std::sqrt((1.0 - P(l, i)*P(l, i))*(1.0 - P(l, k)*P(l, k))) + P(l, i)*P(l, k);
      S(k, i) = p;
      S(i, k) = p;
    }
  }
  return S;
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
// ABC kernels as distributions of observed data y given simulated s (dev/plan.md
// 3.21 (2)): scale is a standard deviation per element; at epsilon = Inf the
// kernel is 1 (log 0), so an SMC can start from the prior predictive
inline arma::colvec ilike_model_abc_scale(const arma::colvec &scale, arma::uword n)
{
  if (scale.n_elem == n) return scale;
  if (scale.n_elem == 1) return arma::colvec(n, arma::fill::value(scale(0)));
  Rcpp::stop("abc kernel: scale has " + std::to_string(scale.n_elem) + " elements, the data " + std::to_string(n) + ".");
}
inline double ilike_model_dabc_uniform(const arma::colvec &y, const arma::colvec &s, double epsilon, const arma::colvec &scale_in, double p)
{
  if (std::isinf(epsilon)) return 0.0;
  if (y.n_elem != s.n_elem)
    Rcpp::stop("abc_uniform: the data have " + std::to_string(y.n_elem) + " elements, the simulated data " + std::to_string(s.n_elem) + ".");
  double n = double(y.n_elem);
  arma::colvec scale = ilike_model_abc_scale(scale_in, y.n_elem);
  arma::colvec d = arma::abs(y - s)/scale;
  double distance = std::isinf(p) ? d.max() : std::pow(arma::accu(arma::pow(d, p)), 1.0/p);
  if (!(distance <= epsilon)) return -std::numeric_limits<double>::infinity();
  // 1 / the volume of the scaled Lp ball of radius epsilon
  double log_unit_ball = std::isinf(p) ? n*std::log(2.0) : n*(std::log(2.0) + std::lgamma(1.0 + 1.0/p)) - std::lgamma(1.0 + n/p);
  return -log_unit_ball - n*std::log(epsilon) - arma::accu(arma::log(scale));
}
inline double ilike_model_dabc_gaussian(const arma::colvec &y, const arma::colvec &s, double epsilon, const arma::colvec &scale_in)
{
  if (std::isinf(epsilon)) return 0.0;
  if (y.n_elem != s.n_elem)
    Rcpp::stop("abc_gaussian: the data have " + std::to_string(y.n_elem) + " elements, the simulated data " + std::to_string(s.n_elem) + ".");
  double n = double(y.n_elem);
  arma::colvec scale = ilike_model_abc_scale(scale_in, y.n_elem);
  arma::colvec z = (y - s)/(epsilon*scale);
  return -0.5*n*std::log(2.0*M_PI) - n*std::log(epsilon) - arma::accu(arma::log(scale)) - 0.5*arma::accu(z % z);
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
      if (this->is_cpp_node(st.variable))
      {
        // x <- cpp::f(a, b, …): the user's arma::mat f(const arma::mat &a, …),
        // each argument a variable (or data, a constant, numbers), read whole
        if ((st.indices.size() != 0) || !st.loops.empty())
          this->error(st.node, "a node computed by C++ (" + st.variable + " <- cpp::f(…)) is a whole variable: no index, no loop.");
        std::string f = st.value["function"].as_string();
        if (find_cpp_function_signature(this->user_code, f).empty())
          this->error(st.value, "cpp::" + f + " is not defined in a cpp section.");
        const RecipeValue &arguments = st.value["arguments"];
        std::string matrix_signature = "inline arma::mat ilike_model_node_matrix_" + st.variable + "(const Parameters &parameters, const Data *data)";
        std::string call = "::" + f + "(";
        std::string reads;
        for (size_t k = 0; k < arguments.size(); ++k)
        {
          if (this->is_deterministic(kind_of(arguments[k]) == "name" ? arguments[k]["name"].as_string() : std::string("")))
            this->error(arguments[k], "the arguments of a C++ node are variables, data, constants or numbers (not other nodes, for now).");
          call += (k ? ", " : "") + this->matrix_code(arguments[k], context);
          reads += this->reads_code(arguments[k], context, true);
        }
        node_code += matrix_signature + ";\n";
        this->deferred_node_code += matrix_signature + "\n{\n  return arma::mat(" + call + "));\n}\n\n";
        std::string reads_signature = "inline void ilike_model_node_reads_" + st.variable + "(std::vector<double> &reads, double term, const Data *data)";
        node_code += reads_signature + ";\n";
        this->deferred_node_code += reads_signature + "\n{\n  Parameters parameters; // (no index read here depends on the parameters)\n" + reads + "}\n\n";
        continue;
      }
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
        Hoist hoist;
        inner_context.hoist = &hoist;
        std::string statement_body = "log_density += " + this->term_code(st, inner_context) + ";\n";
        CodeContext fresh;
        body += "{\n" + hoist.declarations + this->loops_code(st, fresh, statement_body, "") + "}\n";
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

    // the joint prior: its log density, and its simulation in dependency order.
    // (dev/plan.md 3.21, step (b)) With "integrate u with is(…)", u's prior is
    // a separate factor of the inner model (part 1), simulated given the outer
    // parameters, which it may read; the rest is the outer prior (part 0).
    std::vector<size_t> outer_prior, inner_prior;
    std::set<std::string> prior_variables;
    for (size_t k = 0; k < this->model_statements.size(); ++k)
      if (this->model_statements[k].op == "~" && !this->model_statements[k].observed)
      {
        (this->integrated_variables.count(this->model_statements[k].variable) ? inner_prior : outer_prior).push_back(k);
        prior_variables.insert(this->model_statements[k].variable);
      }
    RecipeValue inner_prior_factor;
    std::set<std::string> inner_prior_reads;
    for (int part = 0; part < 2; ++part)
    {
      const bool inner = (part == 1);
      const std::vector<size_t> &prior = inner ? inner_prior : outer_prior;
      const std::string suffix = inner ? "_integrated" : "";
      if (prior.empty())
        continue;
      std::set<std::string> part_variables;
      for (size_t k : prior)
        part_variables.insert(this->model_statements[k].variable);
      // the parameters each statement reads (arguments, indices and loop bounds)
      std::map<size_t, std::set<std::string>> reads;
      bool reads_outer = false; // (the inner prior reads outer parameters)
      for (size_t k : prior)
      {
        const ModelStatement &st = this->model_statements[k];
        std::set<std::string> found;
        this->read_parameters(st.value, found);
        for (size_t i = 0; i < st.indices.size(); ++i) this->read_parameters(st.indices[i], found);
        for (const ModelLoop &loop : st.loops) { this->read_parameters(loop.from, found); this->read_parameters(loop.to, found); }
        for (const ModelLoop &loop : st.loops) found.erase(loop.index);
        for (const std::string &name : found)
        {
          if (!prior_variables.count(name))
            this->error(st.node, "the prior of " + st.variable + " depends on " + name +
                        ", which has no distribution in the model statements (and is not data or a constant).");
          if (!inner && this->integrated_variables.count(name))
            this->error(st.node, "the prior of " + st.variable + " depends on " + name + ", which is integrated out "
                        "(integrate " + name + " with …): only data may depend on a variable integrated out.");
          if (!part_variables.count(name))
            reads_outer = true;
        }
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
      // (the inner prior given the outer parameters starts from them, and
      // returns only its own variables)
      std::string simulate_body = (inner && reads_outer) ? "Parameters parameters = ilike_given;\n" : "Parameters parameters;\n";
      std::set<std::string> allocated;
      for (size_t k : order)
      {
        const ModelStatement &st = this->model_statements[k];
        CodeContext loop_context;
        loop_context.prior = true;
        for (const ModelLoop &loop : st.loops) loop_context.indices.insert(loop.index);
        Hoist hoist;
        CodeContext log_context = loop_context;
        log_context.hoist = &hoist;
        std::string x = this->target_code(st, log_context, false);
        CodeContext fresh;
        fresh.prior = true;
        std::string log_term = "log_density += " + this->density_code(st, x, log_context) + ";\n";
        log_body += "{\n" + hoist.declarations + this->loops_code(st, fresh, log_term, "") + "}\n";

        if (!allocated.count(st.variable))
        {
          allocated.insert(st.variable);
          // the variable's size: the largest index of its statements
          if (!whole_target(st.value["function"].as_string()))
          {
            std::string rows, cols;
            this->variable_extent(st.variable, rows, cols);
            simulate_body += "parameters[" + quoted(st.variable) + "] = arma::mat(" + rows + ", " + cols + ", arma::fill::zeros);\n";
          }
        }
        std::string assign;
        if (whole_target(st.value["function"].as_string()))
          assign = "parameters[" + quoted(st.variable) + "] = " + this->random_code(st, loop_context) + ";\n";
        else
          assign = this->target_code(st, loop_context, true) + " = " + this->random_code(st, loop_context) + ";\n";
        CodeContext fresh2;
        fresh2.prior = true;
        simulate_body += this->loops_code(st, fresh2, assign, "");
      }
      log_body += "return log_density;\n";
      if (inner && reads_outer)
      {
        simulate_body += "Parameters ilike_own;\n";
        for (const std::string &v : part_variables)
          simulate_body += "ilike_own[" + quoted(v) + "] = parameters[" + quoted(v) + "];\n";
        simulate_body += "return ilike_own;\n";
      }
      else
        simulate_body += "return parameters;\n";
      auto indent = [](const std::string &body) {
        std::string out, line;
        std::stringstream lines(body);
        while (std::getline(lines, line)) out += "  " + line + "\n";
        return out;
      };
      this->deferred_function_code += "double ilike_model_log_prior" + suffix + "(const Parameters &parameters)\n{\n  try\n  {\n" + indent(log_body) +
                                      "  }\n" + std::string("  catch (const ilike_model_outside &)\n  {\n    return -std::numeric_limits<double>::infinity(); // (an index that reads parameters is outside its variable)\n  }\n") + "}\n\n";
      RecipeValue factor = RecipeValue::list();
      factor.push_back(this->placeholder("ilike_model_log_prior" + suffix, "evaluate_log_prior"), "evaluate_log_prior");
      if (inner)
      {
        // simulated given the outer parameters (proposal_parameters merged
        // with them, CustomIndependentProposalKernel)
        this->deferred_function_code += "Parameters ilike_model_simulate_prior" + suffix +
                                        "(RandomNumberGenerator &rng, const Parameters &ilike_given)\n{\n" + indent(simulate_body) + "}\n\n";
        factor.push_back(this->placeholder("ilike_model_simulate_prior" + suffix, "simulate_prior_given"), "simulate_prior_given");
      }
      else
      {
        this->deferred_function_code += "Parameters ilike_model_simulate_prior(RandomNumberGenerator &rng)\n{\n" + indent(simulate_body) + "}\n\n";
        factor.push_back(this->placeholder("ilike_model_simulate_prior", "simulate_prior"), "simulate_prior");
      }

      std::vector<const ModelStatement*> statements;
      for (size_t k : order)
        statements.push_back(&this->model_statements[k]);
      if (inner)
      {
        // (no gradient or plate: the inner importance sampler only simulates
        // from this factor, and cancels it)
        inner_prior_factor = factor;
        inner_prior_reads = this->statements_read(statements);
        continue;
      }
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
        std::string many = plate.gradient, many_value;
        if (many.size() > 6 && many.compare(many.size() - 6, 6, "_plate") == 0)
        {
          many_value = many.substr(0, many.size() - 6) + "_many_value_plate";
          many = many.substr(0, many.size() - 6) + "_many_plate";
        }
        else
        {
          many_value = many + "_many_value";
          many = many + "_many";
        }
        entry.push_back(this->placeholder(many, "plate_gradients"), "gradients");
        // and the value with them, in one pass (C222)
        entry.push_back(this->placeholder(many_value, "plate_value_gradients"), "value_gradients");
      }
      RecipeValue variables = RecipeValue::list();
      for (const std::string &name : this->plate_variable_names)
        variables.push_back(RecipeValue(name));
      entry.push_back(variables, "variables");
      this->factor_entries[plate.position].push_back(entry, "plate");
    }

    if (!this->integrations.empty())
      this->integrate_out(inner_prior_factor, inner_prior_reads, outer_prior.empty());

    // a deterministic node defined but never used is harmless; a parameter
    // read by a likelihood but given no distribution may be given by another factor
    this->model_code = code;
  }

  // (dev/plan.md 3.21, step (b)) integrate u with is(points = N): u's prior
  // and the factors of the data that depend on u (directly or through nodes)
  // become an inner model. Its likelihood, ∫ p(y | θ, u) φ(u | θ) du, is
  // estimated for each value of the outer parameters θ by an importance
  // sampler over u whose proposal is φ (which cancels, so the weights are the
  // data's densities): with an ABC kernel as the data's distribution this is
  // the usual ABC likelihood estimate. The inner factors leave the outer list
  // (which is renumbered), and the factor named u, at the integrate
  // statement's place, holds the importance sampler.
  void integrate_out(const RecipeValue &inner_prior_factor, const std::set<std::string> &inner_prior_reads, bool no_outer_prior)
  {
    const Integration &integration = this->integrations[0];
    const std::string &u = integration.variable;
    if (inner_prior_factor.is_null())
      this->error(integration.node, "integrate " + u + ": " + u + " has no distribution in the model statements (" + u + " ~ …).");
    // the data factors that read u, in their order
    std::vector<std::string> inner_observed;
    for (const ModelStatement &st : this->model_statements)
    {
      if (!this->likelihood_term(st) ||
          std::find(inner_observed.begin(), inner_observed.end(), st.variable) != inner_observed.end())
        continue;
      if (this->statements_read({&st}).count(u))
        inner_observed.push_back(st.variable);
    }
    if (inner_observed.empty())
      this->error(integration.node, "integrate " + u + ": no data depend on " + u + ", so there is nothing to integrate.");
    std::sort(inner_observed.begin(), inner_observed.end(), [this](const std::string &a, const std::string &b) {
      return this->model_factor_position.at(a) < this->model_factor_position.at(b);
    });
    RecipeValue inner_factors = RecipeValue::list();
    inner_factors.push_back(inner_prior_factor);
    std::set<size_t> removed;
    // what the factor reads (so that a move on other parameters skips it)
    std::set<std::string> reads = inner_prior_reads;
    for (const std::string &name : inner_observed)
    {
      size_t position = this->model_factor_position.at(name);
      inner_factors.push_back(this->factor_entries[position]);
      removed.insert(position);
      auto found = this->generated_factor_reads.find(position + 1);
      if (found != this->generated_factor_reads.end())
        reads.insert(found->second.begin(), found->second.end());
    }
    reads.erase(u);
    // (a model whose only prior was u's: its outer prior factor is empty)
    if (no_outer_prior && this->model_factor_position.count("model_prior"))
      removed.insert(this->model_factor_position.at("model_prior"));

    // the importance sampler's settings: named (abc(…) is the same as is(…))
    const RecipeValue &arguments = integration.method["arguments"];
    std::vector<std::string> keys = arguments.names();
    RecipeValue parameters = RecipeValue::list();
    for (size_t k = 0; k < arguments.size(); ++k)
    {
      std::string key = (k < keys.size()) ? keys[k] : "";
      if (key.empty())
        this->error(arguments[k], integration.method["function"].as_string() + "(…): its settings are named (e.g. points = 1000).");
      parameters.push_back(this->value_of(arguments[k]), key);
    }
    RecipeValue inner = RecipeValue::list();
    inner.push_back(inner_factors, "factor");
    RecipeValue node = RecipeValue::list();
    node.push_back(RecipeValue("is"), "type");
    node.push_back(RecipeValue(u), "model_file");
    node.push_back(inner, "model");
    node.push_back(parameters, "parameters");
    RecipeValue factor = RecipeValue::list();
    factor.push_back(node, "algorithmic_likelihood");
    size_t position = this->factor_numbers.at(u) - 1;
    this->factor_entries[position] = factor;
    this->generated_factor_reads[position + 1] = reads;
    this->remove_factor_positions(removed);
  }

  // remove factor entries (positions from 0), renumbering the factors after them
  void remove_factor_positions(const std::set<size_t> &removed)
  {
    if (removed.empty())
      return;
    std::vector<size_t> new_position(this->factor_entries.size(), 0);
    std::vector<RecipeValue> kept;
    for (size_t p = 0; p < this->factor_entries.size(); ++p)
      if (!removed.count(p))
      {
        new_position[p] = kept.size();
        kept.push_back(this->factor_entries[p]);
      }
    this->factor_entries = kept;
    for (auto it = this->factor_numbers.begin(); it != this->factor_numbers.end();)
    {
      size_t p = it->second - 1;
      if (removed.count(p))
        it = this->factor_numbers.erase(it);
      else
      {
        it->second = new_position[p] + 1;
        ++it;
      }
    }
    for (auto it = this->model_factor_position.begin(); it != this->model_factor_position.end();)
    {
      if (removed.count(it->second))
        it = this->model_factor_position.erase(it);
      else
      {
        it->second = new_position[it->second];
        ++it;
      }
    }
    std::map<size_t, std::set<std::string>> reads;
    for (const auto &kv : this->generated_factor_reads)
      if (!removed.count(kv.first - 1))
        reads[new_position[kv.first - 1] + 1] = kv.second;
    this->generated_factor_reads = reads;
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
    if (type == "is")
      this->error(call, "is(…) as a factor integrates out a variable of the model statements: "
                  "integrate u with is(points = N).");
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
    if (kind == "sum")
      return this->affine_in(expression["body"], x, why); // (C233: a sum of affine terms; its range reads no parameters)
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
          return {&st, x + "'s distribution is " + f + " (gibbs on a continuous variable needs norm or mvnorm, for a normal full conditional, or gamma, invgamma or beta, for a conjugate one)"};
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

  // Whether an expression is c x^power with c not depending on x (C215):
  // x itself (or an element of it, indices not depending on x) has power 1;
  // products and quotients with x-free expressions keep the power, negation
  // too; 1/e negates it, sqrt(e) halves it, pow(e, k) and e^k with k a number
  // known when compiling multiply it; through deterministic nodes. `why` says
  // why not.
  bool power_of(const RecipeValue &expression, const std::string &x, double &power, std::string &why)
  {
    if (!this->reads_name(expression, x))
    {
      power = 0.0;
      return true;
    }
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
        return this->power_of(this->deterministic_statement(name).value, x, power, why);
      power = 1.0;
      return true;
    }
    if (kind == "unary")
      return this->power_of(expression["operand"], x, power, why);
    double k;
    if (kind == "binary")
    {
      std::string op = expression["operator"].as_string();
      const RecipeValue &left = expression["left"], &right = expression["right"];
      bool l = this->reads_name(left, x), r = this->reads_name(right, x);
      if ((op == "*") || (op == "/"))
      {
        if (l && r)
        {
          why = "both sides of a " + std::string(op == "*" ? "product" : "quotient") + " depend on " + x;
          return false;
        }
        if (!this->power_of(l ? left : right, x, power, why))
          return false;
        if ((op == "/") && r)
          power = -power;
        return true;
      }
      if ((op == "^") && !r && this->compile_time_number(right, k))
      {
        if (!this->power_of(left, x, power, why))
          return false;
        power = power*k;
        return true;
      }
      why = (op == "^") ? "the power of " + x + " is not a number known when compiling"
                        : x + " is in a sum or difference";
      return false;
    }
    if (kind == "call")
    {
      std::string f = expression["function"].as_string();
      const RecipeValue &arguments = expression["arguments"];
      if ((f == "sqrt") && (arguments.size() == 1) && this->power_of(arguments[0], x, power, why))
      {
        power = power/2.0;
        return true;
      }
      if ((f == "pow") && (arguments.size() == 2) && !this->reads_name(arguments[1], x) &&
          this->compile_time_number(arguments[1], k) && this->power_of(arguments[0], x, power, why))
      {
        power = power*k;
        return true;
      }
      why = x + " is inside " + f + "(…)";
      return false;
    }
    why = "it is not a multiple of a power of " + x;
    return false;
  }

  // The conjugate families of gibbs on a positive or (0, 1) variable (C215):
  // x's own distribution gamma, invgamma or beta (arguments not depending on
  // x), and every other term depending on x one whose log density in x is in
  // the family's span, so that log p(x | rest) = a phi1(x) + b phi2(x) + c
  // exactly:
  //   gamma (log x, x): norm/lnorm(m, s) with s = c x^(-1/2) (a precision);
  //     pois(c x); exp(c x); gamma(k, c x) (a rate);
  //   invgamma (log x, 1/x): norm/lnorm(m, s) with s = c x^(1/2) (a
  //     variance); pois(c/x); exp(c/x); gamma(k, c/x);
  //   beta (log x, log(1 - x)): bern(x), binom(n, x).
  // Targets, loops, indices and the other arguments must not depend on x.
  // `family` is gamma, invgamma or beta (or empty: x's distribution is none of
  // these); the answer "" or why not.
  std::pair<const ModelStatement*, std::string> conjugate_problem(const std::string &x, std::string &family)
  {
    family.clear();
    for (const ModelStatement &st : this->model_statements)
      if ((st.op == "~") && !st.observed && (st.variable == x))
      {
        std::string f = st.value["function"].as_string();
        if ((f == "gamma") || (f == "invgamma") || (f == "beta"))
          family = f;
        break;
      }
    if (family.empty())
      return {NULL, ""};
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
        for (size_t k = 0; k < arguments.size(); ++k)
          if (this->reads_name(arguments[k], x))
            return {&st, "the arguments of " + x + "'s distribution depend on " + x};
        continue;
      }
      // the argument through which x enters, and the power it must have
      int which = -1;
      double wanted = 0.0;
      std::string role;
      if (family == "beta")
      {
        if (f == "bern")
          which = 0;
        else if (f == "binom")
          which = 1;
        else
          return {&st, "this term's distribution is " + f + " (with a beta distribution, " + x +
                       " is conjugate for bern(" + x + ") and binom(n, " + x + ") terms)"};
        if (kind_of(arguments[which]) != "name" || arguments[which]["name"].as_string() != x)
        {
          const RecipeValue &a = arguments[which];
          bool element = (kind_of(a) == "index") && (kind_of(a["target"]) == "name") && (a["target"]["name"].as_string() == x);
          if (!element)
            return {&st, "the probability must be " + x + " itself for a beta full conditional"};
        }
      }
      else
      {
        bool gamma = (family == "gamma");
        if ((f == "norm") || (f == "lnorm"))
        {
          which = 1;
          wanted = gamma ? -0.5 : 0.5;
          role = gamma ? "the standard deviation must be c/sqrt(" + x + ") (" + x + " a precision)"
                       : "the standard deviation must be c*sqrt(" + x + ") (" + x + " a variance)";
          if (this->reads_name(arguments[0], x))
            return {&st, x + " is in the mean"};
        }
        else if ((f == "pois") || (f == "exp"))
        {
          which = 0;
          wanted = gamma ? 1.0 : -1.0;
          role = gamma ? "the rate must be c*" + x : "the rate must be c/" + x;
        }
        else if (f == "gamma")
        {
          which = 1;
          wanted = gamma ? 1.0 : -1.0;
          role = gamma ? "the rate must be c*" + x : "the rate must be c/" + x;
          if (this->reads_name(arguments[0], x))
            return {&st, x + " is in the shape"};
        }
        else
          return {&st, "this term's distribution is " + f + " (with a" + std::string(gamma ? " gamma" : "n invgamma") +
                       " distribution, " + x + " is conjugate for norm and lnorm terms through the standard deviation, "
                       "and pois, exp and gamma terms through the rate)"};
        double power;
        std::string why;
        if (!this->power_of(arguments[which], x, power, why))
          return {&st, role + ": " + why};
        if (power != wanted)
          return {&st, role};
      }
    }
    return {NULL, ""};
  }

  // (C228) the Dirichlet full conditional of a probability vector x: x's own
  // distribution dirichlet (arguments free of x), and every other term
  // reading x categorical(x) or multinom(n, x), with x itself as the
  // probabilities and the target, loops and indices free of x
  std::pair<const ModelStatement*, std::string> dirichlet_problem(const std::string &x)
  {
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
        for (size_t k = 0; k < arguments.size(); ++k)
          if (this->reads_name(arguments[k], x))
            return {&st, "the arguments of " + x + "'s distribution depend on " + x};
        continue;
      }
      int which = (f == "categorical") ? 0 : ((f == "multinom") ? 1 : -1);
      if (which < 0)
        return {&st, "this term's distribution is " + f + " (with a dirichlet distribution, " + x +
                     " is conjugate for categorical(" + x + ") and multinom(n, " + x + ") terms)"};
      const RecipeValue &probabilities = arguments[which];
      if ((kind_of(probabilities) != "name") || (probabilities["name"].as_string() != x))
        return {&st, "the probabilities must be " + x + " itself for a Dirichlet full conditional"};
      if ((f == "multinom") && this->reads_name(arguments[0], x))
        return {&st, "the number of trials depends on " + x};
    }
    return {NULL, ""};
  }

  // (C228) the normal-inverse-gamma full conditional of a mean mu and a
  // variance v: v ~ invgamma (arguments free of both), mu ~ norm(m0, c
  // sqrt(v)), and every other term reading either norm(a + b mu, c' sqrt(v))
  // (means affine in mu and free of v; standard deviations c x^(1/2) in v
  // and free of mu), targets, loops and indices free of both
  std::pair<const ModelStatement*, std::string> normal_invgamma_problem(const std::string &mu, const std::string &v)
  {
    bool own_mu = false, own_v = false;
    for (const ModelStatement &st : this->model_statements)
    {
      if (st.op == "<-")
        continue;
      std::set<std::string> reads = this->statements_read(std::vector<const ModelStatement*>(1, &st));
      if (!reads.count(mu) && !reads.count(v))
        continue;
      if (st.op == "+=")
        return {&st, "the term of factor " + st.variable + " depends on " + mu + " or " + v};
      std::string f = st.value["function"].as_string();
      const RecipeValue &arguments = st.value["arguments"];
      for (const ModelLoop &loop : st.loops)
        for (const std::string &x : {mu, v})
          if (this->reads_name(loop.from, x) || this->reads_name(loop.to, x))
            return {&st, "a loop's range depends on " + x};
      for (size_t k = 0; k < st.indices.size(); ++k)
        for (const std::string &x : {mu, v})
          if (this->reads_name(st.indices[k], x))
            return {&st, "an index depends on " + x};
      if (!st.observed && (st.variable == v))
      {
        own_v = true;
        if (f != "invgamma")
          return {&st, v + "'s distribution is " + f + " (the variance needs invgamma)"};
        for (size_t k = 0; k < arguments.size(); ++k)
          if (this->reads_name(arguments[k], mu) || this->reads_name(arguments[k], v))
            return {&st, "the arguments of " + v + "'s distribution depend on " + mu + " or " + v};
        continue;
      }
      if (f != "norm")
        return {&st, "this term's distribution is " + f + " (a normal-inverse-gamma full conditional needs norm terms)"};
      if (!st.observed && (st.variable == mu))
      {
        own_mu = true;
        if (this->reads_name(arguments[0], mu) || this->reads_name(arguments[0], v))
          return {&st, "the mean of " + mu + "'s distribution depends on " + mu + " or " + v};
      }
      else
      {
        std::string why;
        if (this->reads_name(arguments[0], v))
          return {&st, v + " is in the mean"};
        if (!this->affine_in(arguments[0], mu, why))
          return {&st, mu + " enters the mean nonlinearly: " + why};
      }
      double power;
      std::string why;
      if (this->reads_name(arguments[1], mu))
        return {&st, mu + " is in the standard deviation"};
      if (!this->power_of(arguments[1], v, power, why) || (power != 0.5))
        return {&st, "the standard deviation must be c*sqrt(" + v + ") (" + v + " a variance)" + (why.empty() ? "" : ": " + why)};
    }
    if (!own_mu || !own_v)
      return {NULL, (own_mu ? v : mu) + " has no distribution in the model statements"};
    return {NULL, ""};
  }

  // every factor reading x is written as model statements (what the others read is unknown)
  void check_factors_known(const RecipeValue &call, const Arguments &arguments, const std::string &x, const std::string &label)
  {
    if (arguments.has("factors"))
      return;
    std::vector<size_t> numbers = this->automatic_factor_index(x);
    size_t total = this->factor_entries.size();
    for (size_t n = 1; n <= total; ++n)
      if ((numbers.empty() || std::find(numbers.begin(), numbers.end(), n) != numbers.end()) &&
          !this->generated_factor_reads.count(n))
        this->error(call, label + ": factor " + std::to_string(n) + " is not written as model statements, so whether the "
                    "full conditional is conjugate is not known; if it does not depend on " + x +
                    ", list the factors the move uses with factors = [...].");
  }

  // a gibbs move lowered as metropolis_hastings(proposal(…), factors, cut_factors)
  Move gibbs_block_move(const RecipeValue &move_statement, const RecipeValue &call, Arguments &arguments,
                        const std::string &proposal, const RecipeValue &target, const RecipeValue *extra, const std::string &extra_name)
  {
    RecipeValue inner = RecipeValue::list();
    inner.push_back(target, "");
    if (extra != NULL)
      inner.push_back(*extra, extra_name);
    RecipeValue outer = RecipeValue::list();
    outer.push_back(call_with(call, proposal, inner), "");
    for (const std::string &key : {"factors", "cut_factors"})
      if (arguments.has(key))
        outer.push_back(arguments.named[key], key);
    this->lowering_gibbs_block = true;
    Move move = this->lower_move(move_statement, move_statement["target"]["name"].as_string(),
                                 call_with(call, "metropolis_hastings", outer));
    this->lowering_gibbs_block = false;
    return move;
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
    if (gibbs && loops.empty() && (kind_of(element) == "vector"))
    {
      // gibbs([mu, v]): a mean and a variance with a normal-inverse-gamma
      // full conditional (C228)
      const RecipeValue &items = element["elements"];
      if ((items.size() != 2) || (kind_of(items[0]) != "name") || (kind_of(items[1]) != "name"))
        this->error(call, "gibbs on a list takes a mean and a variance: gibbs([mu, v]).");
      std::string a = items[0]["name"].as_string(), b = items[1]["name"].as_string();
      std::string fa, fb;
      for (const ModelStatement &st : this->model_statements)
        if ((st.op == "~") && !st.observed)
        {
          if (st.variable == a) fa = st.value["function"].as_string();
          if (st.variable == b) fb = st.value["function"].as_string();
        }
      bool mean_first = (fb == "invgamma");
      std::string mu = mean_first ? a : b, v = mean_first ? b : a;
      if (!mean_first && (fa != "invgamma"))
        this->error(call, "gibbs([mu, v]): one of the two needs an invgamma distribution (the variance).");
      std::pair<const ModelStatement*, std::string> problem = this->normal_invgamma_problem(mu, v);
      if (!problem.second.empty())
      {
        std::string message = "gibbs([" + a + ", " + b + "]): the full conditional is not normal-inverse-gamma: " + problem.second +
                              ". Use separate gibbs or metropolis moves.";
        if (problem.first != NULL)
          this->error(problem.first->node, message);
        this->error(call, message);
      }
      this->check_factors_known(call, arguments, mu, "gibbs([" + a + ", " + b + "])");
      this->check_factors_known(call, arguments, v, "gibbs([" + a + ", " + b + "])");
      RecipeValue flag = RecipeValue::list();
      flag.push_back(RecipeValue("logical"), "expression");
      flag.push_back(RecipeValue(mean_first), "value");
      return this->gibbs_block_move(move_statement, call, arguments, "normal_invgamma_block", element, &flag, "mean_first");
    }
    if (gibbs && (kind_of(element) != "name") && (kind_of(element) != "index"))
      this->error(call, "gibbs takes a scalar parameter (move m = gibbs(mu)) or an element in a loop "
                  "(for (j in 1:J) move m = gibbs(theta[j])); a proposal such as norm_rw(…) goes in metropolis(…).");
    bool whole_scalar = gibbs && loops.empty() && (kind_of(element) == "name");
    if (whole_scalar)
    {
      // gibbs(mu) on a whole variable: a scalar, as its element mu[1] (C206)
      std::string name = element["name"].as_string();
      bool vector = false, discrete = false;
      for (const ModelStatement &st : this->model_statements)
        if ((st.op == "~") && !st.observed && (st.variable == name))
        {
          std::string f = st.value["function"].as_string();
          if ((st.indices.size() > 0) || whole_target(f))
            vector = true;
          if ((f == "bern") || (f == "binom") || (f == "pois") || (f == "negbinom") || (f == "geom") || (f == "categorical") ||
              (f == "multinom"))
            discrete = true;
        }
      bool dirichlet = false;
      for (const ModelStatement &st : this->model_statements)
        if ((st.op == "~") && !st.observed && (st.variable == name) && (st.value["function"].as_string() == "dirichlet"))
          dirichlet = true;
      if (dirichlet && support.is_null())
      {
        // gibbs(p) on a Dirichlet-distributed probability vector (C228)
        std::pair<const ModelStatement*, std::string> problem = this->dirichlet_problem(name);
        if (!problem.second.empty())
        {
          std::string message = "gibbs(" + name + "): the full conditional of " + name + " is not Dirichlet: " + problem.second +
                                ". Use independent_mh, or a move that keeps it on the simplex.";
          if (problem.first != NULL)
            this->error(problem.first->node, message);
          this->error(call, message);
        }
        this->check_factors_known(call, arguments, name, "gibbs(" + name + ")");
        return this->gibbs_block_move(move_statement, call, arguments, "dirichlet_block", element, NULL, "");
      }
      if (vector)
      {
        // gibbs(beta) on a whole vector (C216): a draw from the block's
        // normal full conditional, if the model statements make it one,
        // lowered as metropolis_hastings(normal_block(beta), …)
        if (discrete || !support.is_null())
          this->error(call, "gibbs(" + name + ") on a whole vector is for a normal full conditional; for the elements of a "
                      "discrete " + name + ", write for (j in 1:J) move m = gibbs(" + name + "[j]).");
        std::pair<const ModelStatement*, std::string> problem = this->normal_conditional_problem(name);
        if (!problem.second.empty())
        {
          std::string message = "gibbs(" + name + "): the full conditional of the vector " + name + " is not normal: " +
                                problem.second + ". Use gibbs on its elements in a loop, if each element's conditional "
                                "is conjugate, or a metropolis move.";
          if (problem.first != NULL)
            this->error(problem.first->node, message);
          this->error(call, message);
        }
        if (!arguments.has("factors"))
        {
          std::vector<size_t> numbers = this->automatic_factor_index(name);
          size_t total = this->factor_entries.size();
          for (size_t n = 1; n <= total; ++n)
            if ((numbers.empty() || std::find(numbers.begin(), numbers.end(), n) != numbers.end()) &&
                !this->generated_factor_reads.count(n))
              this->error(call, "gibbs(" + name + "): factor " + std::to_string(n) + " is not written as model "
                          "statements, so whether its full conditional is normal is not known; if it does not "
                          "depend on " + name + ", list the factors the move uses with factors = [...].");
        }
        RecipeValue inner = RecipeValue::list();
        inner.push_back(element, "");
        RecipeValue proposal_call = call_with(call, "normal_block", inner);
        RecipeValue outer = RecipeValue::list();
        outer.push_back(proposal_call, "");
        for (const std::string &key : {"factors", "cut_factors"})
          if (arguments.has(key))
            outer.push_back(arguments.named[key], key);
        this->lowering_gibbs_block = true;
        Move move = this->lower_move(move_statement, move_statement["target"]["name"].as_string(),
                                     call_with(call, "metropolis_hastings", outer));
        this->lowering_gibbs_block = false;
        return move;
      }
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
            else if (f == "categorical")
            {
              // (C227) 1, ..., K: K from the probabilities, or from the
              // concentrations of a dirichlet-distributed probability vector
              const RecipeValue &probabilities = st.value["arguments"][0];
              size_t K = 0;
              if (this->numeric(substitute(probabilities, constants), number))
                K = number.n_elem;
              else if (kind_of(probabilities) == "name")
                for (const ModelStatement &other : this->model_statements)
                  if ((other.op == "~") && !other.observed && (other.variable == probabilities["name"].as_string()) &&
                      (other.value["function"].as_string() == "dirichlet") &&
                      this->numeric(substitute(other.value["arguments"][0], constants), number))
                    K = number.n_elem;
              for (size_t v = 1; v <= K; ++v)
                states.push_back(double(v));
            }
            break;
          }
        if (states.empty())
        {
          // not discrete: a normal full conditional (C206), or a gamma,
          // inverse gamma or beta one (C215), if the model statements make it one
          std::string family;
          std::pair<const ModelStatement*, std::string> problem = this->conjugate_problem(variable, family);
          if (family.empty())
            problem = this->normal_conditional_problem(variable);
          if (!problem.second.empty())
          {
            std::string message = "gibbs(" + variable + (whole_scalar ? "" : "[…]") + "): " + variable + " is not discrete (bern, binom, or values = [...]), "
                                  "and its full conditional is not " + (family.empty() ? std::string("normal") : family) + ": " + problem.second +
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
                            "statements, so whether its full conditional is conjugate is not known; if it does not "
                            "depend on " + variable + ", list the factors the move uses with factors = [...].");
          }
          type = family.empty() ? "gibbs_normal" : "gibbs_" + family;
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
        if (kind_of(adaptor) != "call" || ((this->builtin_name(adaptor) != "gaussian") && (this->builtin_name(adaptor) != "esjd")))
          this->error(adaptor, "moves on elements are set from SMC particles with smc_adapt = gaussian(c) or esjd(scale = [lower, upper]).");
        if ((this->builtin_name(adaptor) == "esjd") && (type == "flip"))
          this->error(adaptor, "esjd tunes the step of norm_rw and unif_rw on elements; flip has none.");
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
      if ((f == "gamma") || (f == "lnorm") || (f == "exp") || (f == "invgamma") || (f == "halfnorm") ||
          (f == "halfcauchy") || (f == "halft") || (f == "chisq") || (f == "weibull"))
        return transform_item(variable, "log", 0.0, 0.0);
      if (f == "beta")
        return transform_item(variable, "logit", 0.0, 1.0);
      // (C232) whole variables
      if (f == "dirichlet")
        return transform_item(variable, "simplex", 0.0, 0.0);
      if ((f == "wishart") || (f == "invwishart"))
        return transform_item(variable, "cov", 0.0, 0.0);
      if (f == "lkj")
        return transform_item(variable, "corr", 0.0, 0.0);
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
          static const std::string forms = "transform is auto, none, or a list of log(x), logit(x), logit(x, a, b), lower(x, a), upper(x, b), simplex(p), cov(S), corr(R)";
          if (kind_of(item) != "call")
            this->error(item, forms + ".");
          std::string f = item["function"].as_string();
          const RecipeValue &arguments = item["arguments"];
          if (((f != "log") && (f != "logit") && (f != "lower") && (f != "upper") && (f != "simplex") && (f != "cov") && (f != "corr")) ||
              (arguments.size() < 1) || (kind_of(arguments[0]) != "name"))
            this->error(item, forms + ".");
          std::string x = arguments[0]["name"].as_string();
          if (std::find(moved.begin(), moved.end(), x) == moved.end())
            this->error(arguments[0], "transform: " + x + " is not moved by this move.");
          if ((f == "simplex") || (f == "cov") || (f == "corr"))
          {
            if (arguments.size() != 1)
              this->error(item, f + "(" + x + ") takes only the variable.");
            list.push_back(transform_item(x, f, 0.0, 0.0));
            continue;
          }
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

  // a call node like `call`, with another function and arguments (a new node:
  // values share their nodes, so the original is not changed)
  static RecipeValue call_with(const RecipeValue &call, const std::string &function, const RecipeValue &arguments)
  {
    RecipeValue output = RecipeValue::list();
    std::vector<std::string> fields = call.names();
    for (size_t i = 0; i < call.size(); ++i)
    {
      if (fields[i] == "arguments")
        output.push_back(arguments, fields[i]);
      else if (fields[i] == "function")
        output.push_back(RecipeValue(function), fields[i]);
      else
        output.push_back(call[i], fields[i]);
    }
    return output;
  }

  bool lowering_nuts = false;
  size_t sum_temporaries = 0;
  bool lowering_gibbs_block = false;

  // the variables a move changes: the first argument of its proposal (a
  // name, an element or a list), through metropolis(…), choose([…]), …
  void moved_variables(const RecipeValue &call, std::set<std::string> &found) const
  {
    if (kind_of(call) != "call")
      return;
    static const std::set<std::string> wrappers = {"metropolis", "metropolis_hastings", "independent_mh", "unadjusted"};
    std::string f = call["function"].as_string();
    const RecipeValue &arguments = call["arguments"];
    if (arguments.size() == 0)
      return;
    std::vector<std::string> keys = arguments.names();
    if (wrappers.count(f) && call["namespace"].as_string().empty())
    {
      for (size_t k = 0; k < arguments.size(); ++k)
        if ((k < keys.size()) && (keys[k] == "proposal"))
          return this->moved_variables(arguments[k], found);
      return this->moved_variables(arguments[0], found);
    }
    if ((f == "choose") && (kind_of(arguments[0]) == "vector"))
    {
      const RecipeValue &choices = arguments[0]["elements"];
      for (size_t k = 0; k < choices.size(); ++k)
        this->moved_variables(choices[k], found);
      return;
    }
    if (!call["namespace"].as_string().empty() && (call["namespace"].as_string() != "ilike"))
      return; // (a proposal written by the user: what it moves is not known)
    std::function<void(const RecipeValue &)> names = [&](const RecipeValue &e) {
      std::string kind = kind_of(e);
      if (kind == "name")
        found.insert(e["name"].as_string());
      else if (kind == "index")
        names(e["target"]);
      else if (kind == "vector")
        for (size_t k = 0; k < e["elements"].size(); ++k)
          names(e["elements"][k]);
    };
    names(arguments[0]);
  }

  Move lower_move(const RecipeValue &statement, const std::string &name, const RecipeValue &call)
  {
    // (dev/plan.md 3.21, step (b)) a variable integrated out is not a
    // parameter of the outer algorithm
    if (!this->integrated_variables.empty())
    {
      std::set<std::string> moved;
      this->moved_variables(call, moved);
      for (const std::string &v : moved)
        if (this->integrated_variables.count(v))
          this->error(statement, "move " + name + " changes " + v + ", which is integrated out (integrate " + v +
                      " with …), so it is not a parameter here. Remove the integrate statement for inference on " + v + " too.");
    }
    // (dev/plan.md 3.20, S3) a move that changes a parameter of an
    // unnormalised distribution's normaliser would need Z at the proposed value
    if (!this->normaliser_parameters.empty())
    {
      std::set<std::string> moved;
      this->moved_variables(call, moved);
      for (const std::string &v : moved)
      {
        auto found = this->normaliser_parameters.find(v);
        if (found != this->normaliser_parameters.end())
          this->error(statement, "move " + name + " changes " + v + ", a parameter of the normalising constant of " + found->second +
                      "'s distribution, which is not known: the acceptance ratio would need it. Moves on " + found->second +
                      " (with " + v + " fixed) are possible; methods for inferring " + v +
                      " (doubly intractable: exchange, SAV, ABC) are planned (dev/plan.md 3.20).");
      }
    }
    if ((kind_of(call) == "call") && (call["function"].as_string() == "nuts") && call["namespace"].as_string().empty())
    {
      // move m = nuts(x, step = …) (C210): a move of its own, lowered as a
      // Metropolis-Hastings move whose proposal is nuts(…) (its Hamiltonian
      // system); the recipe makes it a NoUTurnMCMC. The move's own arguments
      // stay with the move.
      static const std::set<std::string> move_level = {"factors", "cut_factors", "adapt", "smc_adapt", "transform"};
      const RecipeValue &arguments = call["arguments"];
      std::vector<std::string> keys = arguments.names();
      RecipeValue inner = RecipeValue::list(), outer = RecipeValue::list();
      for (size_t i = 0; i < arguments.size(); ++i)
      {
        std::string key = (i < keys.size()) ? keys[i] : "";
        if (move_level.count(key))
          continue;
        inner.push_back(arguments[i], key);
      }
      outer.push_back(call_with(call, "nuts", inner), "");
      for (size_t i = 0; i < arguments.size(); ++i)
      {
        std::string key = (i < keys.size()) ? keys[i] : "";
        if (move_level.count(key))
          outer.push_back(arguments[i], key);
      }
      this->lowering_nuts = true;
      Move move = this->lower_move(statement, name, call_with(call, "metropolis_hastings", outer));
      this->lowering_nuts = false;
      return move;
    }
    static const std::map<std::string, std::pair<std::string, int>> kinds = {
      {"metropolis", {"m", 3}}, {"metropolis_hastings", {"mh", 1}},
      {"independent_mh", {"imh", 2}}, {"unadjusted", {"unadjusted", 4}}};
    if (kind_of(call) != "call" || !kinds.count(call["function"].as_string()))
      this->error(statement["value"], "a move is metropolis(…), metropolis_hastings(…), independent_mh(…), unadjusted(…), "
                  "nuts(…) or gibbs(…).");
    std::string kind = call["function"].as_string();
    Move move;
    move.name = name;
    move.prefix = kinds.at(kind).first;
    move.code = kinds.at(kind).second;
    const std::string &p = move.prefix;
    // an unadjusted move is biased by any error in the gradient, so it refuses
    // gradients by finite differences (dev/plan.md 3.10; the model section is
    // lowered before the algorithm sections)
    if ((p == "unadjusted") && !this->finite_difference_gradients.empty())
      this->error(call, "unadjusted moves need exact gradients, but the gradient of " + this->finite_difference_gradients[0] +
                        " is computed by finite differences (gradient = autodiff); write the gradient, or use metropolis_hastings.");

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
    if ((kind_of(proposal) == "call") && (proposal["function"].as_string() == "choose") && proposal["namespace"].as_string().empty())
    {
      // choose([q_1, ..., q_K], cost = [...], change = 0.1, smc_adapt = [...])
      // (dev/plan.md 3.19 stage S4, C224): each particle in SMC moves with one
      // of the proposals, the choice tuned from the particles' jumps
      if (p != "mh")
        this->error(proposal, "choose([...]) goes in metropolis_hastings(…).");
      if (arguments.has("smc_adapt"))
        this->error(arguments.named["smc_adapt"], "with choose, give each proposal's smc_adapt inside it: "
                    "choose([...], smc_adapt = [esjd(…), none, …]).");
      if (arguments.has("adapt"))
        this->error(arguments.named["adapt"], "choose's proposals are tuned from the SMC particles (smc_adapt inside choose), not along the chain.");
      Arguments choice = this->split_arguments(proposal);
      this->allow_only(proposal, choice, {"cost", "change", "smc_adapt"}, 1);
      if (choice.positional.empty() || (kind_of(choice.positional[0]) != "vector"))
        this->error(proposal, "choose takes a list of proposals: choose([norm_rw(theta, 1), hmc(theta, step = 0.1, steps = 10)]).");
      const RecipeValue &items = choice.positional[0]["elements"];
      if (items.size() < 2)
        this->error(proposal, "choose needs at least two proposals.");
      RecipeValue components = RecipeValue::list();
      for (size_t k = 0; k < items.size(); ++k)
      {
        const RecipeValue &item = items[k];
        if ((kind_of(item) != "call") || this->is_function_value(item))
          this->error(item, "choose's proposals are built-in proposals (norm_rw, langevin, hmc, …).");
        std::string type = this->builtin_name(item);
        if ((type == "choose") || (type == "nuts") || (type == "normal_block"))
          this->error(item, "choose's proposals are built-in proposals of metropolis_hastings (norm_rw, langevin, hmc, …).");
        if (!role_exists(type, BuiltinRole::MHProposal))
          find_builtin(type, BuiltinRole::MHProposal);
        RecipeValue node = this->builtin_node(item, true, "", true);
        std::string moved = node["variables"].as_string();
        if (variable.empty())
          variable = moved;
        else if (moved != variable)
          this->error(item, "choose's proposals move the same variables.");
        components.push_back(node, "");
      }
      RecipeValue node = RecipeValue::list();
      node.push_back(RecipeValue("choose"), "type");
      node.push_back(RecipeValue(variable), "variables");
      node.push_back(components, "components");
      if (choice.has("smc_adapt"))
      {
        const RecipeValue &list = choice.named["smc_adapt"];
        if ((kind_of(list) != "vector") || (list["elements"].size() != items.size()))
          this->error(list, "choose's smc_adapt lists one adaptor per proposal (none for a proposal without one).");
        RecipeValue adaptors = RecipeValue::list();
        for (size_t k = 0; k < list["elements"].size(); ++k)
        {
          const RecipeValue &item = list["elements"][k];
          if ((kind_of(item) == "name") && (item["name"].as_string() == "none"))
            adaptors.push_back(RecipeValue::list(), "");
          else
            adaptors.push_back(this->builtin_node(item, true, variable), "");
        }
        node.push_back(adaptors, "component_smc_adaptors");
      }
      if (choice.has("cost"))
      {
        arma::mat cost;
        if (!this->numeric(choice.named["cost"], cost) || (cost.n_elem != items.size()) || !(cost.min() > 0.0))
          this->error(choice.named["cost"], "choose's cost lists a positive number per proposal.");
        node.push_back(RecipeValue(std::vector<double>(cost.begin(), cost.end())), "cost");
      }
      if (choice.has("change"))
      {
        double change;
        if (!this->compile_time_number(choice.named["change"], change) || !(change >= 0.0) || !(change <= 1.0))
          this->error(choice.named["change"], "choose's change is a probability.");
        node.push_back(RecipeValue(change), "change");
      }
      entry.push_back(node, p + "_proposal");
      RecipeValue transform = this->move_transform(arguments.has("transform") ? &arguments.named["transform"] : NULL,
                                                   variable, false);
      if (!transform.is_null())
        entry.push_back(transform, "transform");
    }
    else if (kind_of(proposal) == "call" && proposal["function"].as_string() != "density" && !this->is_function_value(proposal))
    {
      std::string type = this->builtin_name(proposal);
      if ((type == "nuts") && !this->lowering_nuts)
        this->error(proposal, "nuts is a move of its own, not a proposal: move m = nuts(x, step = …).");
      if (((type == "normal_block") || (type == "dirichlet_block") || (type == "normal_invgamma_block")) && !this->lowering_gibbs_block)
        this->error(proposal, type + " is not a proposal to write: write move m = gibbs(…).");
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
                                                     variable, (type == "hmc") || (type == "nuts"));
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
      if (this->lowering_nuts)
      {
        // (pretune chooses a number of leapfrog steps, which nuts does not have, C223)
        std::vector<RecipeValue> items;
        if (kind_of(adaptor) == "vector")
          for (size_t k = 0; k < adaptor["elements"].size(); ++k)
            items.push_back(adaptor["elements"][k]);
        else
          items.push_back(adaptor);
        for (const RecipeValue &item : items)
          if ((kind_of(item) == "call") && (this->builtin_name(item) == "pretune"))
            this->error(item, "pretune tunes hmc's step and number of steps; for nuts, use smc_adapt = esjd(step = [lower, upper]).");
      }
      if (kind_of(adaptor) == "vector")
      {
        // several, applied in order: smc_adapt = [gaussian(2.38), esjd(…)] (C211)
        RecipeValue list = RecipeValue::list();
        for (size_t k = 0; k < adaptor["elements"].size(); ++k)
          list.push_back(this->builtin_node(adaptor["elements"][k], true, variable));
        entry.push_back(list, "smc_adaptor");
      }
      else
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
    if (until.is_null() && (kind_of(adaptor) == "call") && (this->builtin_name(adaptor) == "dual_averaging"))
      this->error(adaptor_with_until, "dual_averaging needs until = criterion (e.g. dual_averaging(until = iterations(1000))): "
                  "its iterates do not settle, so the step is fixed at their average when adaptation stops; "
                  "to adapt throughout the run, use scale().");
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
