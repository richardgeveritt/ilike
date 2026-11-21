#include "recipe_internal.h"

// Start-of-run notes (dev/plan.md 3.11): each printed once per run, as a
// message (suppressMessages() hides them) — inner algorithms repeated them
// for every outer point. reset_informs() starts a new run (the R wrappers
// call it).
static std::set<std::string>& informed_notes()
{
  static std::set<std::string> notes;
  return notes;
}

void inform_once(const std::string &note)
{
  if (informed_notes().insert(note).second)
  {
    Rcpp::Function message("message");
    message(note);
  }
}

// [[Rcpp::export]]
void reset_informs()
{
  informed_notes().clear();
}
// Case insensitive string comparison
// From
// https://stackoverflow.com/questions/11635/case-insensitive-string-comparison-in-c

bool ichar_equals(char a, char b) {
  return std::tolower(static_cast<unsigned char>(a)) ==
         std::tolower(static_cast<unsigned char>(b));
}

bool iequals(const std::string &a, const std::string &b) {
  return std::equal(a.begin(), a.end(), b.begin(), b.end(), ichar_equals);
}

std::vector<std::string> split(const std::string &str, const char delimiter) {
  std::vector<std::string> tokens;
  std::stringstream ss(str);
  std::string token;
  while (std::getline(ss, token, delimiter)) {
    tokens.push_back(token);
  }
  return tokens;
}

bool isDoubleInList(const RecipeValue &list, const std::string &name) {
  RecipeValue var = list[name];
  return (var).is_real() && (var).size() == 1;
}

bool isNumericVectorInList(const RecipeValue &list, const std::string &name) {
  RecipeValue var = list[name];
  return (var).is_vector() && (var).is_numeric();
}

bool isNumericMatrixInList(const RecipeValue &list, const std::string &name) {
  RecipeValue var = list[name];
  return (var).is_matrix() && (var).is_numeric();
}

Data get_data(const RecipeValue &model) {
  if (model.contains("data")) {
    RecipeValue data_list = model["data"];
    RecipeValue data_list0 = data_list[0];
    RecipeValue data_value = data_list0["data"];

    Data output = (data_value).function_pointer<DataPtr>()();
    return output;

  } else {
    return Data();
  }
}

// The arguments of the built-in simulation-based likelihoods and filters,
// bound by the registry (builtins.cpp; dev/plan.md 3.14, I2c) and returned
// in the order the factor-building code reads them.

RecipeValue get_abc_euclidean_uniform_parameter_info(const RecipeValue &model_parameters,
                                              const RecipeValue &current_sbi,
                                              const std::string &sbi_name) {
  BuiltinArguments a = bind_builtin(current_sbi, BuiltinRole::SBILikelihood, model_parameters);
  std::string scale_variable = a.has("scale_variable") ? a.variable("scale_variable") : std::string("");
  arma::colvec scale = a.has("scale") ? a.vector("scale") : arma::colvec();
  return RecipeValue::list({int(a.integer("number_of_points")), a.variable("tolerance_variable"),
                            a.real("tolerance"), a.logical("store_output"), scale_variable, scale,
                            a.logical("parallel"), int(a.integer("grain_size"))});
}

RecipeValue get_abc_lp_uniform_parameter_info(const RecipeValue &model_parameters,
                                       const RecipeValue &current_sbi,
                                       const std::string &sbi_name) {
  BuiltinArguments a = bind_builtin(current_sbi, BuiltinRole::SBILikelihood, model_parameters);
  std::string scale_variable = a.has("scale_variable") ? a.variable("scale_variable") : std::string("");
  arma::colvec scale = a.has("scale") ? a.vector("scale") : arma::colvec();
  return RecipeValue::list({int(a.integer("number_of_points")), a.variable("tolerance_variable"),
                            a.real("tolerance"), a.real("p"), a.logical("store_output"),
                            scale_variable, scale, a.logical("parallel"),
                            int(a.integer("grain_size"))});
}

RecipeValue get_abc_enki_parameter_info(const RecipeValue &model_parameters,
                                 const RecipeValue &current_sbi,
                                 const std::string &sbi_name) {
  BuiltinArguments a = bind_builtin(current_sbi, BuiltinRole::SBILikelihood, model_parameters);
  std::string scale_variable = a.has("scale_variable") ? a.variable("scale_variable") : std::string("");
  arma::colvec scale = a.has("scale") ? a.vector("scale") : arma::colvec();
  return RecipeValue::list({int(a.integer("number_of_points")), int(a.integer("estimator_type")),
                            a.variable("tolerance_variable"), a.real("epsilon"), a.string("shifter"),
                            size_t(a.integer("number_of_targets")), a.real("significance_level"),
                            int(a.integer("lag")), scale_variable, scale, a.logical("on_summary"),
                            a.logical("parallel"), int(a.integer("grain_size"))});
}

RecipeValue get_sl_parameter_info(const RecipeValue &model_parameters,
                           const RecipeValue &current_sbi,
                           const std::string &sbi_name) {
  BuiltinArguments a = bind_builtin(current_sbi, BuiltinRole::SBILikelihood, model_parameters);
  return RecipeValue::list({int(a.integer("number_of_points")), a.logical("unbiased"),
                            a.logical("store_output"), a.logical("parallel"),
                            int(a.integer("grain_size"))});
}

// an algorithmic likelihood's call, with the inner model's file (kept by
// the parser as model_file) as its first argument
static BuiltinArguments bind_algorithmic_likelihood(const RecipeValue &current_likelihood,
                                                    const RecipeValue &model_parameters) {
  RecipeValue node = RecipeValue::list();
  node.push_back(current_likelihood["type"], "type");
  if (current_likelihood.contains("model_file"))
    node.push_back(current_likelihood["model_file"], "variables");
  if (current_likelihood.contains("parameters"))
    node.push_back(current_likelihood["parameters"], "parameters");
  return bind_builtin(node, BuiltinRole::AlgorithmicLikelihood, model_parameters);
}

RecipeValue get_kf_info(const RecipeValue &model_parameters, const RecipeValue &current_likelihood,
                 const std::string &likelihood_name) {
  BuiltinArguments a = bind_algorithmic_likelihood(current_likelihood, model_parameters);
  return RecipeValue::list({int(a.integer("iterations_to_store"))});
}

RecipeValue get_enkf_info(const RecipeValue &model_parameters, const RecipeValue &current_likelihood,
                   const std::string &likelihood_name) {
  BuiltinArguments a = bind_algorithmic_likelihood(current_likelihood, model_parameters);
  return RecipeValue::list({size_t(a.integer("ensemble_members")), a.string("shifter"),
                            int(a.integer("iterations_to_store")), a.logical("parallel"),
                            size_t(a.integer("grain_size"))});
}

// integrate u with smc(particles = N, …) (dev/plan.md 3.21, (c1))
RecipeValue get_integrate_smc_info(const RecipeValue &model_parameters, const RecipeValue &current_likelihood) {
  BuiltinArguments a = bind_algorithmic_likelihood(current_likelihood, model_parameters);
  return RecipeValue::list({size_t(a.integer("particles")), a.logical("parallel"), size_t(a.integer("grain_size"))});
}

// integrate u with is(points = N) (dev/plan.md 3.21, step (b))
RecipeValue get_is_info(const RecipeValue &model_parameters, const RecipeValue &current_likelihood) {
  BuiltinArguments a = bind_algorithmic_likelihood(current_likelihood, model_parameters);
  return RecipeValue::list({size_t(a.integer("points")), a.logical("parallel"), size_t(a.integer("grain_size"))});
}

RecipeValue get_pf_info(const RecipeValue &model_parameters, const RecipeValue &current_likelihood,
                 const std::string &likelihood_name) {
  BuiltinArguments a = bind_algorithmic_likelihood(current_likelihood, model_parameters);
  // iterations_to_store defaults to 1: keep the final generation (its points
  // and weights) for each outer point (owner, 2026-11-15); grain_size to 64,
  // as PF() (dev/plan.md 2.6b, G5)
  return RecipeValue::list({size_t(a.integer("particles")), int(a.integer("iterations_to_store")),
                            a.logical("parallel"), size_t(a.integer("grain_size"))});
}

std::vector<size_t> get_enk_likelihood_indices(const RecipeValue &model,
                                               const RecipeValue &model_parameters) {
  std::vector<size_t> enk_likelihood_indices;
  if (model.contains("method")) {
    RecipeValue method_info = model["method"];

    for (size_t i = 0; i < method_info.size(); ++i) {
      RecipeValue current_method = method_info[i];
      if (current_method.contains("enk_likelihood_index")) {
        std::vector<double> indices;

        RecipeValue index_info_value = current_method["enk_likelihood_index"];
        indices = Rcpp::as<std::vector<double>>((index_info_value).function_pointer<NumericVectorPtr>()());

        enk_likelihood_indices.reserve(indices.size());
        for (size_t j = 0; j < indices.size(); ++j) {
          enk_likelihood_indices.push_back(indices[j]-1);
        }

        return enk_likelihood_indices;
      }
    }
  } else {
    return enk_likelihood_indices;
  }

  return enk_likelihood_indices;
}

RecipeValue get_ssm_info(const RecipeValue &model, const RecipeValue &model_parameters) {
  if (!model.contains("method"))
    stop("Method not found in model.");
  RecipeValue method_info = model["method"];
  for (size_t i = 0; i < method_info.size(); ++i) {
    RecipeValue current_method = method_info[i];
    if (current_method.contains("ssm")) {
      BuiltinArguments a = bind_builtin(current_method["ssm"], BuiltinRole::SSM, model_parameters);
      return RecipeValue::list({a.variable("index"), size_t(a.integer("first_index")),
                                size_t(a.integer("last_index")), a.variable("time"),
                                a.real("initial_time"), a.variable("time_difference"),
                                a.real("time_step"), size_t(a.integer("predictions_per_update"))});
    }
  }
  stop("SSM not found in model.");
}

bool is_linear_gaussian_data_model(const RecipeValue &current_factor) {
  return current_factor.contains("linear_gaussian_data_covariance") ||
         current_factor.contains("linear_gaussian_data_model");
}

bool is_linear_gaussian_transition_model(const RecipeValue &current_factor) {
  return current_factor.contains(
             "linear_gaussian_transition_covariance") ||
         current_factor.contains("linear_gaussian_transition_model");
}

// The one-line built-in: (first variable, [second variable,] matrix, covariance).
RecipeValue get_one_line_linear_gaussian(const RecipeValue &current_factor,
                                  const std::string &block_name,
                                  size_t number_of_variables,
                                  const RecipeValue &model_parameters) {
  if (!(current_factor[block_name]).is_list())
    stop(block_name + ": invalid specification.");
  BuiltinRole role = (number_of_variables == 2) ? BuiltinRole::LinearGaussianDataModel
                                                : BuiltinRole::LinearGaussianTransitionModel;
  BuiltinArguments a = bind_builtin(current_factor[block_name], role, model_parameters);
  std::string first_variable = a.variable("variable");
  std::string second_variable =
      (number_of_variables == 2) ? a.string("state_variables") : first_variable;
  return RecipeValue::list({first_variable, second_variable, a.matrix("matrix"),
                            a.matrix("covariance")});
}

LinearGaussianSpec get_linear_gaussian_data_model(const RecipeValue &current_factor,
                                                  const RecipeValue &model_parameters) {
  LinearGaussianSpec spec;
  if (current_factor.contains("linear_gaussian_data_model")) {
    RecipeValue info = get_one_line_linear_gaussian(
        current_factor, "linear_gaussian_data_model", 2, model_parameters);
    spec.variable = (info[0]).as_string();
    spec.state_variable_string = (info[1]).as_string();
    spec.matrix = (info[2]).as_matrix();
    spec.covariance = (info[3]).as_matrix();
  } else {
    if (!current_factor.contains("linear_gaussian_data_matrix"))
      stop("linear_gaussian_data must contain both matrix and covariance.");
    if (!current_factor.contains("type"))
      stop("linear_gaussian_data must contain a type.");
    if (!current_factor.contains("linear_gaussian_data_variable"))
      stop("linear_gaussian_data must contain linear_gaussian_data_variable.");
    if (!current_factor.contains(
            "linear_gaussian_data_state_variable"))
      stop("linear_gaussian_data must contain "
           "linear_gaussian_data_state_variable.");
    RecipeValue matrix_value = current_factor["linear_gaussian_data_matrix"];
    RecipeValue covariance_value = current_factor["linear_gaussian_data_covariance"];
    spec.variable = (current_factor["linear_gaussian_data_variable"]).function_pointer<StringPtr>()();
    spec.state_variable_string =
        (current_factor["linear_gaussian_data_state_variable"]).function_pointer<StringPtr>()();
    size_t type = current_factor["type"];
    if (type == 1) {
      spec.matrix = (matrix_value).function_pointer<GetNoParamsMatrixPtr>()();
      spec.covariance = (covariance_value).function_pointer<GetNoParamsMatrixPtr>()();
    } else if (type == 2) {
      spec.is_function = true;
      spec.matrix_function = (matrix_value).function_pointer<GetMatrixPtr>();
      spec.covariance_function = (covariance_value).function_pointer<GetMatrixPtr>();
    } else {
      stop("linear_gaussian_data specificiation is invalid.");
    }
  }
  spec.state_variables = split(spec.state_variable_string, ';');
  return spec;
}

LinearGaussianSpec
get_linear_gaussian_transition_model(const RecipeValue &current_factor,
                                     const RecipeValue &model_parameters,
                                     bool variable_required) {
  LinearGaussianSpec spec;
  if (current_factor.contains("linear_gaussian_transition_model")) {
    RecipeValue info = get_one_line_linear_gaussian(
        current_factor, "linear_gaussian_transition_model", 1,
        model_parameters);
    spec.variable = (info[0]).as_string();
    spec.matrix = (info[2]).as_matrix();
    spec.covariance = (info[3]).as_matrix();
  } else {
    if (!current_factor.contains("linear_gaussian_transition_matrix"))
      stop("linear_gaussian must contain both matrix and covariance.");
    if (!current_factor.contains("type"))
      stop("linear_gaussian must contain a type.");
    if (current_factor.contains(
            "linear_gaussian_transition_variable"))
      spec.variable =
          (current_factor["linear_gaussian_transition_variable"]).function_pointer<StringPtr>()();
    else if (variable_required)
      stop("linear_gaussian must contain linear_gaussian_transition_variable.");
    RecipeValue matrix_value = current_factor["linear_gaussian_transition_matrix"];
    RecipeValue covariance_value =
        current_factor["linear_gaussian_transition_covariance"];
    size_t type = current_factor["type"];
    if (type == 1) {
      spec.matrix = (matrix_value).function_pointer<GetNoParamsMatrixPtr>()();
      spec.covariance = (covariance_value).function_pointer<GetNoParamsMatrixPtr>()();
    } else if (type == 2) {
      spec.is_function = true;
      spec.matrix_function = (matrix_value).function_pointer<GetMatrixPtr>();
      spec.covariance_function = (covariance_value).function_pointer<GetMatrixPtr>();
    } else {
      stop("linear_gaussian specificiation is invalid.");
    }
  }
  spec.state_variable_string = spec.variable;
  spec.state_variables = std::vector<std::string>(1, spec.variable);
  return spec;
}
