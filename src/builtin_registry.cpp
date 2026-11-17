// The registry of built-ins (dev/plan.md 3.14, I2): binding the arguments of
// a call to its entry, and construction. The entries are in builtins.cpp.

#include <RcppArmadillo.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <sstream>

#include "builtin_registry.h"
#include "symmetric_proposal_kernel.h"

namespace ilike
{

// --- roles and types ---------------------------------------------------------

std::string role_block(BuiltinRole role)
{
  switch (role)
  {
    case BuiltinRole::Prior: return "prior";
    case BuiltinRole::ISProposal: return "is_proposal";
    case BuiltinRole::IMHProposal: return "imh_proposal";
    case BuiltinRole::MHProposal: return "mh_proposal";
    case BuiltinRole::UnadjustedProposal: return "unadjusted_proposal";
    case BuiltinRole::MProposal: return "m_proposal";
    case BuiltinRole::TransitionModel: return "transition_model";
    case BuiltinRole::TransitionProposal: return "transition_proposal";
    case BuiltinRole::MCMCAdaptor: return "mcmc_adaptor";
    case BuiltinRole::SMCAdaptor: return "smc_adaptor";
    case BuiltinRole::MCMCTermination: return "mcmc_termination";
    case BuiltinRole::AdaptiveResampling: return "adaptive_resampling";
    case BuiltinRole::AdaptiveTarget: return "adaptive_target";
    case BuiltinRole::SMCTermination: return "smc_termination";
    case BuiltinRole::EnKShifter: return "enk_shifter";
    case BuiltinRole::SMCSequence: return "smc_sequence";
    case BuiltinRole::EnKSequence: return "enk_sequence";
    case BuiltinRole::SBILikelihood: return "sbi_likelihood";
    case BuiltinRole::AlgorithmicLikelihood: return "algorithmic_likelihood";
    case BuiltinRole::SSM: return "ssm";
    case BuiltinRole::LinearGaussianDataModel: return "linear_gaussian_data_model";
    case BuiltinRole::LinearGaussianTransitionModel: return "linear_gaussian_transition_model";
  }
  return "";
}

std::string role_description(BuiltinRole role)
{
  switch (role)
  {
    case BuiltinRole::Prior: return "prior";
    case BuiltinRole::ISProposal: return "importance sampling proposal";
    case BuiltinRole::IMHProposal: return "independent MH proposal";
    case BuiltinRole::MHProposal: return "MH proposal";
    case BuiltinRole::UnadjustedProposal: return "unadjusted proposal";
    case BuiltinRole::MProposal: return "Metropolis proposal";
    case BuiltinRole::TransitionModel: return "transition model";
    case BuiltinRole::TransitionProposal: return "transition proposal";
    case BuiltinRole::MCMCAdaptor: return "MCMC adaptor";
    case BuiltinRole::SMCAdaptor: return "SMC adaptor";
    case BuiltinRole::MCMCTermination: return "MCMC termination";
    case BuiltinRole::AdaptiveResampling: return "adaptive resampling method";
    case BuiltinRole::AdaptiveTarget: return "adaptive target method";
    case BuiltinRole::SMCTermination: return "SMC termination";
    case BuiltinRole::EnKShifter: return "EnK shifter";
    case BuiltinRole::SMCSequence: return "SMC sequence";
    case BuiltinRole::EnKSequence: return "EnK sequence";
    case BuiltinRole::SBILikelihood: return "simulation-based likelihood";
    case BuiltinRole::AlgorithmicLikelihood: return "algorithmic likelihood";
    case BuiltinRole::SSM: return "state-space model setting";
    case BuiltinRole::LinearGaussianDataModel: return "linear-Gaussian data model";
    case BuiltinRole::LinearGaussianTransitionModel: return "linear-Gaussian transition model";
  }
  return "";
}

std::string argument_type_name(ArgumentType type)
{
  switch (type)
  {
    case ArgumentType::Variable: return "variable";
    case ArgumentType::Real: return "real";
    case ArgumentType::Integer: return "integer";
    case ArgumentType::Logical: return "logical";
    case ArgumentType::String: return "string";
    case ArgumentType::Vector: return "vector";
    case ArgumentType::Matrix: return "matrix";
    case ArgumentType::Flag: return "flag";
  }
  return "";
}

ArgumentSpec &ArgumentSpec::between(double lower_in, double upper_in, bool open)
{
  this->lower = lower_in;
  this->upper = upper_in;
  this->lower_open = open;
  this->upper_open = open;
  return *this;
}

ArgumentSpec &ArgumentSpec::at_least(double lower_in)
{
  this->lower = lower_in;
  this->lower_open = false;
  return *this;
}

ArgumentSpec &ArgumentSpec::above(double lower_in)
{
  this->lower = lower_in;
  this->lower_open = true;
  return *this;
}

ArgumentSpec required(const std::string &name, ArgumentType type, const std::string &description)
{
  ArgumentSpec spec;
  spec.name = name;
  spec.type = type;
  spec.required = true;
  spec.description = description;
  return spec;
}

ArgumentSpec optional(const std::string &name, ArgumentType type, const RecipeValue &default_value,
                      const std::string &description)
{
  ArgumentSpec spec = required(name, type, description);
  spec.required = false;
  spec.default_value = default_value;
  return spec;
}

ArgumentSpec optional(const std::string &name, ArgumentType type, const std::string &description)
{
  return optional(name, type, RecipeValue(), description);
}

ArgumentSpec flag(const std::string &name, const std::string &description)
{
  return optional(name, ArgumentType::Flag, RecipeValue(false), description);
}

bool Builtin::has_role(BuiltinRole role) const
{
  return std::find(this->roles.begin(), this->roles.end(), role) != this->roles.end();
}

// --- bound arguments ---------------------------------------------------------

BuiltinArguments::BuiltinArguments(const std::string &builtin_in, BuiltinRole role_in)
  : builtin_name(builtin_in), builtin_role(role_in)
{
}

std::string BuiltinArguments::where() const
{
  return "ilike::" + this->builtin_name + " (" + role_block(this->builtin_role) + ")";
}

bool BuiltinArguments::given(const std::string &name) const
{
  return this->given_names.count(name) > 0;
}

bool BuiltinArguments::has(const std::string &name) const
{
  return this->values.count(name) > 0;
}

void BuiltinArguments::set(const std::string &name, const RecipeValue &value, bool given_in)
{
  this->values[name] = value;
  if (given_in)
    this->given_names.insert(name);
}

const RecipeValue &BuiltinArguments::get(const std::string &name) const
{
  auto found = this->values.find(name);
  if (found == this->values.end())
    Rcpp::stop(this->where() + ": the argument " + name + " has no value.");
  return found->second;
}

std::string BuiltinArguments::variable(const std::string &name) const { return this->get(name).as_string(); }
std::string BuiltinArguments::string(const std::string &name) const { return this->get(name).as_string(); }
double BuiltinArguments::real(const std::string &name) const { return this->get(name).as_double(); }
long BuiltinArguments::integer(const std::string &name) const { return long(this->get(name).as_double()); }
bool BuiltinArguments::logical(const std::string &name) const { return this->get(name).as_bool(); }
arma::colvec BuiltinArguments::vector(const std::string &name) const { return this->get(name).as_colvec(); }
arma::mat BuiltinArguments::matrix(const std::string &name) const { return this->get(name).as_matrix(); }

// --- reading arguments given as text -----------------------------------------

static std::string trim(const std::string &text)
{
  size_t first = text.find_first_not_of(" \t\r\n");
  if (first == std::string::npos)
    return "";
  size_t last = text.find_last_not_of(" \t\r\n");
  return text.substr(first, last - first + 1);
}

static bool is_identifier(const std::string &text)
{
  if (text.empty() || !(std::isalpha((unsigned char)text[0]) || text[0] == '_' || text[0] == '.'))
    return false;
  for (char c : text)
    if (!(std::isalnum((unsigned char)c) || c == '_' || c == '.'))
      return false;
  return true;
}

// the position of a top-level "=" that is not part of "==", "<=", ">=", "!="
static size_t top_level_assignment(const std::string &text)
{
  int depth = 0;
  char quote = 0;
  for (size_t i = 0; i < text.size(); ++i)
  {
    char c = text[i];
    if (quote)
    {
      if (c == '\\') { ++i; continue; }
      if (c == quote) quote = 0;
      continue;
    }
    if (c == '"' || c == '\'') { quote = c; continue; }
    if (c == '(' || c == '[' || c == '{') ++depth;
    else if (c == ')' || c == ']' || c == '}') --depth;
    else if (c == '=' && depth == 0)
    {
      bool before = (i > 0) && (text[i-1] == '=' || text[i-1] == '<' || text[i-1] == '>' || text[i-1] == '!');
      bool after = (i + 1 < text.size()) && (text[i+1] == '=');
      if (!before && !after)
        return i;
    }
  }
  return std::string::npos;
}

static bool parse_number(const std::string &text, double &value)
{
  if (text.empty())
    return false;
  const char *start = text.c_str();
  char *end = nullptr;
  value = std::strtod(start, &end);
  return (end != start) && (*end == '\0');
}

static bool parse_logical(const std::string &text, bool &value)
{
  if (text == "TRUE" || text == "T" || text == "true") { value = true; return true; }
  if (text == "FALSE" || text == "F" || text == "false") { value = false; return true; }
  return false;
}

static bool is_quoted(const std::string &text)
{
  return (text.size() >= 2) &&
         ((text.front() == '"' && text.back() == '"') || (text.front() == '\'' && text.back() == '\''));
}

// p1, p2, ...: the k-th run-time constant (old files; L-6)
static bool positional_constant(const std::string &text, size_t &index)
{
  if (text.size() < 2 || text[0] != 'p')
    return false;
  for (size_t i = 1; i < text.size(); ++i)
    if (!std::isdigit((unsigned char)text[i]))
      return false;
  index = std::stoul(text.substr(1));
  return index >= 1;
}

static std::string ordinal(size_t k)
{
  std::string suffix = "th";
  if ((k % 100) < 11 || (k % 100) > 13)
  {
    if (k % 10 == 1) suffix = "st";
    else if (k % 10 == 2) suffix = "nd";
    else if (k % 10 == 3) suffix = "rd";
  }
  return std::to_string(k) + suffix;
}

static std::string given_count(const RecipeValue &constants)
{
  size_t n = constants.size();
  return std::to_string(n) + (n == 1 ? " was given" : " were given");
}

static std::string constant_names(const RecipeValue &constants)
{
  std::string names;
  std::vector<std::string> keys = constants.names();
  for (const std::string &key : keys)
    if (!key.empty())
      names = names + (names.empty() ? "" : ", ") + key;
  return names.empty() ? "none by name" : names;
}

// --- canonical calls ---------------------------------------------------------

// an argument given as text ("2", "scale=2") is split into its name and
// value; one given as a value (from the ilike language) keeps the name it has
static void push_argument(RecipeValue &arguments, const RecipeValue &raw, const std::string &key = "")
{
  if (!key.empty())
  {
    arguments.push_back(raw.is_string() && raw.size() == 1 ? RecipeValue(trim(raw.as_string())) : raw, key);
    return;
  }
  if (raw.is_string() && raw.size() == 1)
  {
    std::string text = trim(raw.as_string());
    // nothing between the parentheses, as in ilike::averaged()
    if (text.empty())
      return;
    size_t equals = top_level_assignment(text);
    if (equals != std::string::npos)
    {
      std::string name = trim(text.substr(0, equals));
      if (is_identifier(name))
      {
        arguments.push_back(RecipeValue(trim(text.substr(equals + 1))), name);
        return;
      }
    }
    arguments.push_back(RecipeValue(text));
    return;
  }
  arguments.push_back(raw);
}

RecipeValue canonical_builtin_call(const RecipeValue &node)
{
  if (!node.is_list())
    Rcpp::stop("Recipe: a built-in call must be a list.");
  if (node.contains("builtin"))
    return node;

  RecipeValue call = RecipeValue::list();
  RecipeValue arguments = RecipeValue::list();
  if (node.contains("type"))
  {
    call.push_back(node["type"], "builtin");
    if (node.contains("variables"))
      push_argument(arguments, node["variables"]);
    if (node.contains("parameters"))
    {
      const RecipeValue &parameters = node["parameters"];
      std::vector<std::string> keys = parameters.names();
      if (parameters.is_list())
        for (size_t i = 0; i < parameters.size(); ++i)
          push_argument(arguments, parameters[i], i < keys.size() ? keys[i] : "");
      else if (!parameters.is_null())
        push_argument(arguments, parameters);
    }
  }
  else if (node.contains("method"))
  {
    call.push_back(node["method"], "builtin");
    if (node.contains("values"))
    {
      const RecipeValue &values = node["values"];
      std::vector<std::string> keys = values.names();
      if (values.is_list())
        for (size_t i = 0; i < values.size(); ++i)
          push_argument(arguments, values[i], i < keys.size() ? keys[i] : "");
      else if (!values.is_null())
        push_argument(arguments, values);
    }
  }
  else
  {
    Rcpp::stop("Recipe: a built-in call needs the name of the built-in.");
  }
  call.push_back(arguments, "arguments");
  return call;
}

// --- lookup ------------------------------------------------------------------

static bool same_name(const std::string &a, const std::string &b)
{
  if (a.size() != b.size())
    return false;
  for (size_t i = 0; i < a.size(); ++i)
    if (std::tolower((unsigned char)a[i]) != std::tolower((unsigned char)b[i]))
      return false;
  return true;
}

const Builtin &find_builtin(const std::string &name, BuiltinRole role)
{
  const std::vector<Builtin> &registry = builtin_registry();
  const Builtin *other_role = nullptr;
  for (const Builtin &entry : registry)
  {
    if (!same_name(entry.name, name))
      continue;
    if (entry.has_role(role))
      return entry;
    other_role = &entry;
  }
  if (other_role != nullptr)
  {
    std::string blocks;
    for (BuiltinRole r : other_role->roles)
      blocks = blocks + (blocks.empty() ? "" : ", ") + role_block(r);
    Rcpp::stop("ilike::" + other_role->name + " cannot be used as a" +
               (std::string("aeiouAEIOU").find(role_description(role)[0]) != std::string::npos ? "n " : " ") +
               role_description(role) + " (" + role_block(role) + "); it can be used in: " + blocks + ".");
  }
  std::string available;
  for (const Builtin &entry : registry)
    if (entry.has_role(role))
      available = available + (available.empty() ? "" : ", ") + "ilike::" + entry.name;
  Rcpp::stop("Unknown " + role_description(role) + " ilike::" + name +
             " (available: " + available + ").");
}

// --- binding -----------------------------------------------------------------

static std::string value_text(const RecipeValue &value)
{
  if (value.is_null())
    return "";
  if (value.kind() == RecipeValue::Kind::Logical)
    return value.as_bool() ? "TRUE" : "FALSE";
  if (value.is_string())
    return "\"" + value.as_string() + "\"";
  if (value.is_numeric() && value.size() == 1)
  {
    std::ostringstream out;
    out << value.as_double();
    return out.str();
  }
  return value.kind_name();
}

static std::string signature(const Builtin &entry)
{
  std::string text;
  for (const ArgumentSpec &spec : entry.arguments)
  {
    std::string item = spec.name;
    if (spec.type == ArgumentType::Flag)
      item = item + " (flag)";
    else if (!spec.required && !spec.default_value.is_null())
      item = item + " = " + value_text(spec.default_value);
    else if (!spec.required)
      item = item + " (optional)";
    text = text + (text.empty() ? "" : ", ") + item;
  }
  return "ilike::" + entry.name + "(" + text + ")";
}

static RecipeValue resolve_argument(const RecipeValue &raw, const ArgumentSpec &spec,
                                    const RecipeValue &constants, const std::string &where)
{
  std::string context = where + ", argument " + spec.name;
  RecipeValue value = raw;

  if (raw.is_string() && raw.size() == 1)
  {
    std::string text = trim(raw.as_string());
    size_t index = 0;

    if (spec.type == ArgumentType::Variable)
    {
      if (is_quoted(text))
        text = text.substr(1, text.size() - 2);
      if (text.empty())
        Rcpp::stop(context + ": a variable name is needed.");
      if (text.find(';') != std::string::npos)
        Rcpp::stop(context + ": only one variable is allowed.");
      return RecipeValue(text);
    }

    if (spec.type == ArgumentType::String)
    {
      if (positional_constant(text, index) && !constants.contains(text))
      {
        if (index > constants.size())
          Rcpp::stop(context + ": " + text + " refers to the " + ordinal(index) +
                     " run-time constant, but " + given_count(constants) + ".");
        value = constants[index - 1];
      }
      else
      {
        if (is_quoted(text))
          text = text.substr(1, text.size() - 2);
        return RecipeValue(text);
      }
    }
    else
    {
      double number;
      bool logical;
      if (parse_number(text, number))
        value = RecipeValue(number);
      else if (parse_logical(text, logical))
        value = RecipeValue::logicals(std::vector<int>(1, logical ? 1 : 0));
      else if (is_identifier(text) && constants.contains(text))
        value = constants[text];
      else if (positional_constant(text, index))
      {
        if (index > constants.size())
          Rcpp::stop(context + ": " + text + " refers to the " + ordinal(index) +
                     " run-time constant, but " + given_count(constants) + ".");
        value = constants[index - 1];
      }
      else if (is_identifier(text))
        Rcpp::stop(context + ": " + text + " is not the name of a run-time constant (constants given: " +
                   constant_names(constants) + ").");
      else
        Rcpp::stop(context + ": \"" + text + "\" is not a " + argument_type_name(spec.type) +
                   " or the name of a run-time constant.");
    }
  }

  // convert to the argument's type
  switch (spec.type)
  {
    case ArgumentType::String:
      if (!value.is_string() || value.size() != 1)
        Rcpp::stop(context + ": a string is needed, found " + value.kind_name() + ".");
      return value;
    case ArgumentType::Variable:
      if (!value.is_string() || value.size() != 1)
        Rcpp::stop(context + ": a variable name is needed, found " + value.kind_name() + ".");
      return value;
    case ArgumentType::Logical:
    case ArgumentType::Flag:
      if (!value.is_numeric() || value.size() != 1)
        Rcpp::stop(context + ": TRUE or FALSE is needed, found " + value.kind_name() + ".");
      return RecipeValue::logicals(std::vector<int>(1, value.as_bool() ? 1 : 0));
    case ArgumentType::Real:
    case ArgumentType::Integer:
    {
      if (!value.is_numeric() || value.size() != 1)
        Rcpp::stop(context + ": a number is needed, found " + value.kind_name() +
                   (value.size() != 1 ? " of length " + std::to_string(value.size()) : "") + ".");
      double number = value.as_double();
      if (spec.type == ArgumentType::Integer && (number != std::floor(number)))
        Rcpp::stop(context + ": a whole number is needed, found " + value_text(RecipeValue(number)) + ".");
      bool below = spec.lower_open ? !(number > spec.lower) : !(number >= spec.lower);
      bool beyond = spec.upper_open ? !(number < spec.upper) : !(number <= spec.upper);
      if (below || beyond)
      {
        std::ostringstream range;
        if (std::isfinite(spec.lower) && std::isfinite(spec.upper))
          range << "between " << spec.lower << " and " << spec.upper
                << (spec.lower_open ? " (exclusive)" : "");
        else if (std::isfinite(spec.lower))
          range << (spec.lower_open ? "above " : "at least ") << spec.lower;
        else
          range << (spec.upper_open ? "below " : "at most ") << spec.upper;
        Rcpp::stop(context + ": must be " + range.str() + ", found " + value_text(RecipeValue(number)) + ".");
      }
      return RecipeValue(number);
    }
    case ArgumentType::Vector:
      if (!value.is_numeric())
        Rcpp::stop(context + ": a numeric vector is needed, found " + value.kind_name() + ".");
      return RecipeValue(value.as_doubles());
    case ArgumentType::Matrix:
      if (value.is_matrix())
        return value;
      if (value.is_numeric() && value.size() == 1)
        return RecipeValue(value.as_matrix());
      Rcpp::stop(context + ": a matrix is needed, found " + value.kind_name() +
                 (value.size() != 1 ? " of length " + std::to_string(value.size()) : "") + ".");
  }
  return value;
}

BuiltinArguments bind_builtin(const RecipeValue &node, BuiltinRole role,
                              const RecipeValue &constants)
{
  RecipeValue call = canonical_builtin_call(node);
  const Builtin &entry = find_builtin(call["builtin"].as_string(), role);
  BuiltinArguments bound(entry.name, role);
  std::string where = bound.where();

  const RecipeValue &arguments = call["arguments"];
  std::vector<std::string> keys = arguments.names();

  auto find_spec = [&](const std::string &name) -> const ArgumentSpec* {
    for (const ArgumentSpec &spec : entry.arguments)
      if (spec.name == name)
        return &spec;
    return nullptr;
  };
  auto bind = [&](const ArgumentSpec &spec, const RecipeValue &raw) {
    if (bound.given(spec.name))
      Rcpp::stop(where + ": the argument " + spec.name + " is given more than once.");
    bound.set(spec.name, resolve_argument(raw, spec, constants, where), true);
  };

  // keyword arguments first, then positional ones fill the remaining
  // arguments in order (as in R); a flag may also be given by its name alone
  std::vector<size_t> positional;
  for (size_t i = 0; i < arguments.size(); ++i)
  {
    const std::string key = (i < keys.size()) ? keys[i] : "";
    if (key.empty())
    {
      const RecipeValue &raw = arguments[i];
      if (raw.is_string() && raw.size() == 1)
      {
        const ArgumentSpec *spec = find_spec(trim(raw.as_string()));
        if ((spec != nullptr) && (spec->type == ArgumentType::Flag))
        {
          bind(*spec, RecipeValue(true));
          continue;
        }
      }
      positional.push_back(i);
      continue;
    }
    const ArgumentSpec *spec = find_spec(key);
    if (spec == nullptr)
      Rcpp::stop(where + ": no argument called " + key + " (it takes " + signature(entry) + ").");
    bind(*spec, arguments[i]);
  }

  size_t next = 0;
  for (size_t i : positional)
  {
    while ((next < entry.arguments.size()) &&
           (bound.given(entry.arguments[next].name) || entry.arguments[next].type == ArgumentType::Flag))
      ++next;
    if (next >= entry.arguments.size())
      Rcpp::stop(where + ": too many arguments (it takes " + signature(entry) + ").");
    bind(entry.arguments[next], arguments[i]);
    ++next;
  }

  for (const ArgumentSpec &spec : entry.arguments)
  {
    if (bound.given(spec.name))
      continue;
    if (spec.required)
      Rcpp::stop(where + ": the argument " + spec.name + " is missing (it takes " + signature(entry) + ").");
    if (!spec.default_value.is_null())
      bound.set(spec.name, resolve_argument(spec.default_value, spec, constants, where), false);
  }
  return bound;
}

// --- construction ------------------------------------------------------------

static void missing_constructor(const Builtin &entry, BuiltinRole role)
{
  Rcpp::stop("Registry: ilike::" + entry.name + " has no constructor for " + role_block(role) + ".");
}

DistributionFactor *make_builtin_prior(const RecipeValue &node, const RecipeValue &constants)
{
  BuiltinArguments arguments = bind_builtin(node, BuiltinRole::Prior, constants);
  const Builtin &entry = find_builtin(arguments.builtin(), BuiltinRole::Prior);
  if (!entry.prior)
    missing_constructor(entry, BuiltinRole::Prior);
  return entry.prior(arguments);
}

IndependentProposalKernel *make_builtin_independent_proposal(const RecipeValue &node, BuiltinRole role,
                                                             const RecipeValue &constants)
{
  BuiltinArguments arguments = bind_builtin(node, role, constants);
  const Builtin &entry = find_builtin(arguments.builtin(), role);
  if (!entry.independent_proposal)
    missing_constructor(entry, role);
  return entry.independent_proposal(arguments);
}

ProposalKernel *make_builtin_proposal(const RecipeValue &node, BuiltinRole role,
                                      const RecipeValue &constants)
{
  BuiltinArguments arguments = bind_builtin(node, role, constants);
  const Builtin &entry = find_builtin(arguments.builtin(), role);
  if (entry.proposal)
    return entry.proposal(arguments);
  if (entry.symmetric_proposal)
    return static_cast<ProposalKernel*>(entry.symmetric_proposal(arguments));
  missing_constructor(entry, role);
  return nullptr;
}

SymmetricProposalKernel *make_builtin_symmetric_proposal(const RecipeValue &node, BuiltinRole role,
                                                         const RecipeValue &constants)
{
  BuiltinArguments arguments = bind_builtin(node, role, constants);
  const Builtin &entry = find_builtin(arguments.builtin(), role);
  if (!entry.symmetric_proposal)
    missing_constructor(entry, role);
  return entry.symmetric_proposal(arguments);
}

MCMCAdaptor *make_builtin_mcmc_adaptor(const RecipeValue &node, const MCMCAdaptorContext &context,
                                       const RecipeValue &constants)
{
  BuiltinArguments arguments = bind_builtin(node, BuiltinRole::MCMCAdaptor, constants);
  const Builtin &entry = find_builtin(arguments.builtin(), BuiltinRole::MCMCAdaptor);
  if (!entry.mcmc_adaptor)
    missing_constructor(entry, BuiltinRole::MCMCAdaptor);
  return entry.mcmc_adaptor(arguments, context);
}

SMCAdaptor *make_builtin_smc_adaptor(const RecipeValue &node, ProposalKernel *proposal,
                                     const RecipeValue &constants)
{
  BuiltinArguments arguments = bind_builtin(node, BuiltinRole::SMCAdaptor, constants);
  const Builtin &entry = find_builtin(arguments.builtin(), BuiltinRole::SMCAdaptor);
  if (!entry.smc_adaptor)
    missing_constructor(entry, BuiltinRole::SMCAdaptor);
  return entry.smc_adaptor(arguments, proposal);
}

MCMCTermination *make_builtin_mcmc_termination(const RecipeValue &node, const RecipeValue &constants)
{
  BuiltinArguments arguments = bind_builtin(node, BuiltinRole::MCMCTermination, constants);
  const Builtin &entry = find_builtin(arguments.builtin(), BuiltinRole::MCMCTermination);
  if (!entry.mcmc_termination)
    missing_constructor(entry, BuiltinRole::MCMCTermination);
  return entry.mcmc_termination(arguments);
}

SMCCriterion *make_builtin_smc_criterion(const RecipeValue &node, BuiltinRole role,
                                         const RecipeValue &constants, size_t number_of_particles)
{
  BuiltinArguments arguments = bind_builtin(node, role, constants);
  const Builtin &entry = find_builtin(arguments.builtin(), role);
  if (!entry.smc_criterion)
    missing_constructor(entry, role);
  return entry.smc_criterion(arguments, number_of_particles);
}

SMCTermination *make_builtin_smc_termination(const RecipeValue &node, const RecipeValue &constants)
{
  BuiltinArguments arguments = bind_builtin(node, BuiltinRole::SMCTermination, constants);
  const Builtin &entry = find_builtin(arguments.builtin(), BuiltinRole::SMCTermination);
  if (!entry.smc_termination)
    missing_constructor(entry, BuiltinRole::SMCTermination);
  return entry.smc_termination(arguments);
}

EnsembleShifter *make_builtin_enk_shifter(const RecipeValue &node, const RecipeValue &constants)
{
  BuiltinArguments arguments = bind_builtin(node, BuiltinRole::EnKShifter, constants);
  const Builtin &entry = find_builtin(arguments.builtin(), BuiltinRole::EnKShifter);
  if (!entry.enk_shifter)
    missing_constructor(entry, BuiltinRole::EnKShifter);
  return entry.enk_shifter(arguments);
}

EnsembleShifter *make_builtin_enk_shifter(const std::string &name)
{
  RecipeValue node = RecipeValue::list();
  node.push_back(RecipeValue(name), "builtin");
  node.push_back(RecipeValue::list(), "arguments");
  return make_builtin_enk_shifter(node, RecipeValue::list());
}

// --- the table, for the documentation and tests -------------------------------

RecipeValue builtin_table()
{
  RecipeValue table = RecipeValue::list();
  for (const Builtin &entry : builtin_registry())
  {
    std::vector<std::string> roles;
    for (BuiltinRole role : entry.roles)
      roles.push_back(role_block(role));
    RecipeValue arguments = RecipeValue::list();
    for (const ArgumentSpec &spec : entry.arguments)
    {
      RecipeValue argument = RecipeValue::list();
      argument.push_back(RecipeValue(spec.name), "name");
      argument.push_back(RecipeValue(argument_type_name(spec.type)), "type");
      argument.push_back(RecipeValue(spec.required), "required");
      argument.push_back(RecipeValue(value_text(spec.default_value)), "default");
      argument.push_back(RecipeValue(spec.description), "description");
      arguments.push_back(argument, spec.name);
    }
    RecipeValue item = RecipeValue::list();
    item.push_back(RecipeValue(entry.name), "name");
    item.push_back(RecipeValue(entry.description), "description");
    item.push_back(RecipeValue(roles), "roles");
    item.push_back(arguments, "arguments");
    table.push_back(item, entry.name);
  }
  return table;
}

}
