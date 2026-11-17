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

// from
// https://stackoverflow.com/questions/2165921/converting-from-a-stdstring-to-bool
bool stob(const std::string &s) {
  auto result = false; // failure to assert is false

  std::istringstream is(s);
  // first try simple integer conversion
  is >> result;

  if (is.fail()) {
    // simple integer failed; try boolean
    is.clear();
    is >> std::boolalpha >> result;
  }

  if (is.fail()) {
    stop(s + " is not convertable to bool");
  }

  return result;
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

bool isDoubleInList(const Rcpp::List &list, int index_from_one) {
  if (index_from_one <= list.size()) {
    SEXP var = list[index_from_one - 1];
    return Rf_isReal(var) && Rf_length(var) == 1;
  } else {
    return false;
  }
}

bool isStringInList(const Rcpp::List &list, int index_from_one) {
  if (index_from_one <= list.size()) {
    SEXP var = list[index_from_one - 1];
    return Rf_isString(var) && Rf_length(var) == 1;
  } else {
    return false;
  }
}

bool isVectorInList(const Rcpp::List &list, int index_from_one) {
  if (index_from_one <= list.size()) {
    SEXP var = list[index_from_one - 1];
    return Rf_isVector(var);
  } else {
    return false;
  }
}

bool isMatrixInList(const Rcpp::List &list, int index_from_one) {
  if (index_from_one <= list.size()) {
    SEXP var = list[index_from_one - 1];
    return Rf_isMatrix(var);
  } else {
    return false;
  }
}

bool isNumericVectorInList(const Rcpp::List &list, int index_from_one) {
  if (index_from_one <= list.size()) {
    SEXP var = list[index_from_one - 1];
    return Rf_isVector(var) && Rf_isNumeric(var);
  } else {
    return false;
  }
}

bool isNumericMatrixInList(const Rcpp::List &list, int index_from_one) {
  if (index_from_one <= list.size()) {
    SEXP var = list[index_from_one - 1];
    return Rf_isMatrix(var) && Rf_isNumeric(var);
  } else {
    return false;
  }
}

bool isStringVectorInList(const Rcpp::List &list, const std::string &name) {
  SEXP var = list[name];
  return Rf_isVector(var) && Rf_isString(var);
}

bool isDoubleInList(const Rcpp::List &list, const std::string &name) {
  SEXP var = list[name];
  return Rf_isReal(var) && Rf_length(var) == 1;
}

bool isNumericVectorInList(const Rcpp::List &list, const std::string &name) {
  SEXP var = list[name];
  return Rf_isVector(var) && Rf_isNumeric(var);
}

bool isNumericMatrixInList(const Rcpp::List &list, const std::string &name) {
  SEXP var = list[name];
  return Rf_isMatrix(var) && Rf_isNumeric(var);
}

Data get_data(const List &model) {
  if (model.containsElementNamed("data")) {
    List data_list = model["data"];
    List data_list0 = data_list[0];
    SEXP data_SEXP = data_list0["data"];

    Data output = load_data(data_SEXP);
    return output;

  } else {
    return Data();
  }
}

bool extract_bool_parameter(const List &parameters_from_file,
                            const List &model_parameters, size_t index) {
  bool result;
  std::string parameter_string = parameters_from_file[index];

  if (parameter_string.at(0) == 'p') {
    size_t parameter_index;
    try {
      parameter_index = std::stoi(parameter_string.substr(1));
    } catch (...) {
      Rcpp::stop("Parameter index in model file not an integer (did you supply "
                 "the appropriate number of model parameters?).");
    }

    if (isDoubleInList(model_parameters, parameter_index)) {
      result = model_parameters[parameter_index - 1];
    } else {
      Rcpp::stop("extract_bool_parameter: Parameter index does not correspond "
                 "to a real number in the parameters file.");
    }
  } else {
    try {
      result = stob(parameter_string);
    } catch (...) {
      Rcpp::stop(
          "Parameter in model file is not a real number - something is wrong "
          "with a list of arguments you provided in the ilike file.");
    }
  }

  return result;
}

int extract_int_parameter(const List &parameters_from_file,
                          const List &model_parameters, size_t index) {
  int result;
  std::string parameter_string = parameters_from_file[index];

  if (parameter_string.at(0) == 'p') {
    size_t parameter_index;
    try {
      parameter_index = std::stoi(parameter_string.substr(1));
    } catch (...) {
      Rcpp::stop("Parameter index in model file not an integer (did you supply "
                 "the appropriate number of model parameters?).");
    }

    if (isDoubleInList(model_parameters, parameter_index)) {
      result = model_parameters[parameter_index - 1];
    } else {
      Rcpp::stop("extract_int_parameter: Parameter index does not correspond "
                 "to a real number in the parameters file. Did you supply a "
                 "model_parameter_list argument?");
    }
  } else {
    try {
      result = std::stoi(parameter_string);
    } catch (...) {
      Rcpp::stop(
          "Parameter in model file is not a real number - something is wrong "
          "with a list of arguments you provided in the ilike file.");
    }
  }

  return result;
}

double extract_double_parameter(const List &parameters_from_file,
                                const List &model_parameters, size_t index) {
  double result;
  std::string parameter_string = parameters_from_file[index];

  if (parameter_string.at(0) == 'p') {
    size_t parameter_index;
    try {
      parameter_index = std::stoi(parameter_string.substr(1));
    } catch (...) {
      Rcpp::stop("Parameter index in model file not an integer (did you supply "
                 "the appropriate number of model parameters?).");
    }

    if (isDoubleInList(model_parameters, parameter_index)) {
      result = model_parameters[parameter_index - 1];
    } else {
      Rcpp::stop("extract_double_parameter: Parameter index does not "
                 "correspond to a real number in the parameters file.");
    }
  } else {
    try {
      result = std::stod(parameter_string);
    } catch (...) {
      Rcpp::stop(
          "Parameter in model file is not a real number - something is wrong "
          "with a list of arguments you provided in the ilike file.");
    }
  }

  return result;
}

std::string extract_string_parameter(const List &parameters_from_file,
                                     const List &model_parameters,
                                     size_t index) {
  std::string result;
  std::string parameter_string = parameters_from_file[index];

  if (parameter_string.at(0) == 'p') {
    size_t parameter_index;
    try {
      parameter_index = std::stoi(parameter_string.substr(1));
    } catch (...) {
      std::string temp_string = parameters_from_file[index];
      result = temp_string;
    }

    if (isStringInList(model_parameters, parameter_index)) {
      std::string temp_string = model_parameters[parameter_index - 1];
      result = temp_string;
    } else {
      Rcpp::stop("extract_string_parameter: Parameter index does not "
                 "correspond to a string in the parameters file.");
    }
  } else {
    try {
      std::string temp_string = parameters_from_file[index];
      result = temp_string;
    } catch (...) {
      Rcpp::stop("Parameter in model file is not a string - something is wrong "
                 "with a list of arguments you provided in the ilike file.");
    }
  }

  return result;
}

arma::colvec extract_vector_parameter(const List &parameters_from_file,
                                      const List &model_parameters,
                                      size_t index) {
  arma::colvec result;
  std::string parameter_string = parameters_from_file[index];

  if (parameter_string.at(0) == 'p') {
    size_t parameter_index;
    try {
      parameter_index = std::stoi(parameter_string.substr(1));
    } catch (...) {
      Rcpp::stop("Parameter index in model file not an integer (did you supply "
                 "the appropriate number of model parameters?).");
    }

    if (isVectorInList(model_parameters, parameter_index)) {
      result = as<NumericVector>(model_parameters[parameter_index - 1]);
    } else {
      Rcpp::stop("extract_vector_parameter: Parameter index does not "
                 "correspond to a real number in the parameters file.");
    }
  } else {
    try {
      result = std::stod(parameter_string);
    } catch (...) {
      Rcpp::stop(
          "Parameter in model file is not a real number - something is wrong "
          "with a list of arguments you provided in the ilike file.");
    }
  }

  return result;
}

arma::mat extract_matrix_parameter(const List &parameters_from_file,
                                   const List &model_parameters, size_t index) {
  arma::mat result;
  std::string parameter_string = parameters_from_file[index];

  if (parameter_string.at(0) == 'p') {
    size_t parameter_index;
    try {
      parameter_index = std::stoi(parameter_string.substr(1));
    } catch (...) {
      Rcpp::stop("Parameter index in model file not an integer (did you supply "
                 "the appropriate number of model parameters?).");
    }

    if (isMatrixInList(model_parameters, parameter_index)) {
      result = as<arma::mat>(model_parameters[parameter_index - 1]);
    } else {
      Rcpp::stop("extract_matrix_parameter: Parameter index does not "
                 "correspond to a real number in the parameters file.");
    }
  } else {
    try {
      result = std::stod(parameter_string);
    } catch (...) {
      Rcpp::stop(
          "Parameter in model file is not a real number - something is wrong "
          "with a list of arguments you provided in the ilike file.");
    }
  }

  return result;
}

List get_single_variable_one_double_parameter_info(
    const List &model_parameters, const List &current_distribution,
    const std::string &distribution_name) {
  std::string augmented_variable_names = current_distribution["variables"];

  std::vector<std::string> variable_names =
      split(augmented_variable_names, ';');

  // throw error if more than one variable
  if (variable_names.size() != 1)
    Rcpp::stop("Only one variable allowed for " + distribution_name + ".");

  if (!current_distribution.containsElementNamed("parameters"))
    Rcpp::stop("Missing parameters for " + distribution_name +
               " (one parameters required).");

  List parameters = current_distribution["parameters"];

  if (parameters.size() != 1)
    Rcpp::stop("One parameter required for " + distribution_name + ".");

  double first_param =
      extract_double_parameter(parameters, model_parameters, 0);

  return List::create(variable_names[0], first_param);
}

List get_single_variable_two_double_parameter_info(
    const List &model_parameters, const List &current_distribution,
    const std::string &distribution_name) {
  std::string augmented_variable_names = current_distribution["variables"];

  std::vector<std::string> variable_names =
      split(augmented_variable_names, ';');

  // throw error if more than one variable
  if (variable_names.size() != 1)
    Rcpp::stop("Only one variable allowed for " + distribution_name + ".");

  if (!current_distribution.containsElementNamed("parameters"))
    Rcpp::stop("Missing parameters for " + distribution_name +
               " (two parameters required).");

  List parameters = current_distribution["parameters"];

  if (parameters.size() != 2)
    Rcpp::stop("Two parameters required for " + distribution_name + ".");

  double first_param =
      extract_double_parameter(parameters, model_parameters, 0);

  double second_param =
      extract_double_parameter(parameters, model_parameters, 1);

  return List::create(variable_names[0], first_param, second_param);
}

List get_single_variable_vector_parameter_info(
    const List &model_parameters, const List &current_distribution,
    const std::string &distribution_name) {
  std::string augmented_variable_names = current_distribution["variables"];

  std::vector<std::string> variable_names =
      split(augmented_variable_names, ';');

  // throw error if more than one variable
  if (variable_names.size() != 1)
    Rcpp::stop("Only one variable allowed for " + distribution_name + ".");

  if (!current_distribution.containsElementNamed("parameters"))
    Rcpp::stop("Missing parameters for " + distribution_name +
               " (one parameter required).");

  List parameters = current_distribution["parameters"];

  if (parameters.size() != 1)
    Rcpp::stop("One parameter required for " + distribution_name + ".");

  arma::colvec first_param =
      extract_vector_parameter(parameters, model_parameters, 0);

  return List::create(variable_names[0], first_param);
}

List get_single_variable_matrix_parameter_info(
    const List &model_parameters, const List &current_distribution,
    const std::string &distribution_name) {
  std::string augmented_variable_names = current_distribution["variables"];

  std::vector<std::string> variable_names =
      split(augmented_variable_names, ';');

  // throw error if more than one variable
  if (variable_names.size() != 1)
    Rcpp::stop("Only one variable allowed for " + distribution_name + ".");

  if (!current_distribution.containsElementNamed("parameters"))
    Rcpp::stop("Missing parameters for " + distribution_name +
               " (one parameter required).");

  List parameters = current_distribution["parameters"];

  if (parameters.size() != 1)
    Rcpp::stop("One parameter required for " + distribution_name + ".");

  arma::mat first_param =
      extract_matrix_parameter(parameters, model_parameters, 0);

  return List::create(variable_names[0], first_param);
}

List get_single_variable_vector_and_matrix_parameter_info(
    const List &model_parameters, const List &current_distribution,
    const std::string &distribution_name) {
  std::string augmented_variable_names = current_distribution["variables"];

  std::vector<std::string> variable_names =
      split(augmented_variable_names, ';');

  // throw error if more than one variable
  if (variable_names.size() != 1)
    Rcpp::stop("Only one variable allowed for " + distribution_name + ".");

  if (!current_distribution.containsElementNamed("parameters"))
    Rcpp::stop("Missing parameters for " + distribution_name +
               " (two parameters required).");

  List parameters = current_distribution["parameters"];

  if (parameters.size() != 2)
    Rcpp::stop("Two parameters required for " + distribution_name + ".");

  arma::colvec first_param =
      extract_vector_parameter(parameters, model_parameters, 0);

  arma::mat second_param =
      extract_matrix_parameter(parameters, model_parameters, 1);

  return List::create(variable_names[0], first_param, second_param);
}

List get_single_variable_matrix_and_double_parameter_info(
    const List &model_parameters, const List &current_distribution,
    const std::string &distribution_name) {
  std::string augmented_variable_names = current_distribution["variables"];

  std::vector<std::string> variable_names =
      split(augmented_variable_names, ';');

  // throw error if more than one variable
  if (variable_names.size() != 1)
    Rcpp::stop("Only one variable allowed for " + distribution_name + ".");

  if (!current_distribution.containsElementNamed("parameters"))
    Rcpp::stop("Missing parameters for " + distribution_name +
               " (two parameters required).");

  List parameters = current_distribution["parameters"];

  if (parameters.size() != 2)
    Rcpp::stop("Two parameters required for " + distribution_name + ".");

  arma::mat first_param =
      extract_matrix_parameter(parameters, model_parameters, 0);

  double second_param =
      extract_double_parameter(parameters, model_parameters, 1);

  return List::create(variable_names[0], first_param, second_param);
}

// (variable, matrix) or (variable, matrix, double): the double is optional
// and takes default_double when it is not given.
List get_single_variable_matrix_and_optional_double_parameter_info(
    const List &model_parameters, const List &current_distribution,
    const std::string &distribution_name, double default_double) {
  std::string augmented_variable_names = current_distribution["variables"];
  std::vector<std::string> variable_names =
      split(augmented_variable_names, ';');

  if (variable_names.size() != 1)
    Rcpp::stop("Only one variable allowed for " + distribution_name + ".");

  if (!current_distribution.containsElementNamed("parameters"))
    Rcpp::stop("Missing parameters for " + distribution_name +
               " (a matrix and optionally a scale required).");

  List parameters = current_distribution["parameters"];

  if ((parameters.size() != 1) && (parameters.size() != 2))
    Rcpp::stop(distribution_name +
               " requires a matrix and optionally a scale (one or two "
               "parameters after the variable).");

  arma::mat first_param =
      extract_matrix_parameter(parameters, model_parameters, 0);

  double second_param = default_double;
  if (parameters.size() == 2)
    second_param = extract_double_parameter(parameters, model_parameters, 1);

  return List::create(variable_names[0], first_param, second_param);
}

List get_single_variable_vector_and_matrix_and_double_parameter_info(
    const List &model_parameters, const List &current_distribution,
    const std::string &distribution_name) {
  std::string augmented_variable_names = current_distribution["variables"];

  std::vector<std::string> variable_names =
      split(augmented_variable_names, ';');

  // throw error if more than one variable
  if (variable_names.size() != 1)
    Rcpp::stop("Only one variable allowed for " + distribution_name + ".");

  if (!current_distribution.containsElementNamed("parameters"))
    Rcpp::stop("Missing parameters for " + distribution_name +
               " (three parameters required).");

  List parameters = current_distribution["parameters"];

  if (parameters.size() != 3)
    Rcpp::stop("Three parameters required for " + distribution_name + ".");

  arma::colvec first_param =
      extract_vector_parameter(parameters, model_parameters, 0);

  arma::mat second_param =
      extract_matrix_parameter(parameters, model_parameters, 1);

  double third_param =
      extract_double_parameter(parameters, model_parameters, 2);

  return List::create(variable_names[0], first_param, second_param,
                      third_param);
}

List get_abc_euclidean_uniform_parameter_info(const List &model_parameters,
                                              const List &current_sbi,
                                              const std::string &sbi_name) {

  if (!current_sbi.containsElementNamed("parameters"))
    Rcpp::stop("Missing parameters for " + sbi_name +
               " (4-8 parameters required).");

  List parameters = current_sbi["parameters"];

  if (!((parameters.size() == 3) || (parameters.size() == 4) ||
        (parameters.size() == 5) || (parameters.size() == 6) ||
        (parameters.size() == 7) || (parameters.size() == 8)))
    Rcpp::stop("3-8 parameters required for " + sbi_name + ".");

  int number_of_points = extract_int_parameter(parameters, model_parameters, 0);

  std::string tolerance_variable =
      extract_string_parameter(parameters, model_parameters, 1);

  double tolerance = extract_double_parameter(parameters, model_parameters, 2);

  bool store_output = false;
  if (parameters.size() > 3) {
    store_output = extract_bool_parameter(parameters, model_parameters, 3);
  }

  std::string scale_variable = "";
  if (parameters.size() > 4) {
    scale_variable = extract_string_parameter(parameters, model_parameters, 4);
  }

  arma::colvec scale = arma::colvec();
  if (parameters.size() > 5) {
    scale = extract_vector_parameter(parameters, model_parameters, 5);
  }

  bool parallel = false;
  if (parameters.size() > 6) {
    parallel = extract_bool_parameter(parameters, model_parameters, 6);
  }

  int grain_size = 100000;
  if (parameters.size() > 7) {
    grain_size = extract_int_parameter(parameters, model_parameters, 7);
  }

  return List::create(number_of_points, tolerance_variable, tolerance,
                      store_output, scale_variable, scale, parallel,
                      grain_size);
}

List get_abc_lp_uniform_parameter_info(const List &model_parameters,
                                       const List &current_sbi,
                                       const std::string &sbi_name) {

  if (!current_sbi.containsElementNamed("parameters"))
    Rcpp::stop("Missing parameters for " + sbi_name +
               " (4-9 parameters required).");

  List parameters = current_sbi["parameters"];

  if (!((parameters.size() == 4) || (parameters.size() == 5) ||
        (parameters.size() == 6) || (parameters.size() == 7) ||
        (parameters.size() == 8) || (parameters.size() == 9)))
    Rcpp::stop("4-9 parameters required for " + sbi_name + ".");

  int number_of_points = extract_int_parameter(parameters, model_parameters, 0);

  std::string tolerance_variable =
      extract_string_parameter(parameters, model_parameters, 1);

  double tolerance = extract_double_parameter(parameters, model_parameters, 2);

  double p = extract_double_parameter(parameters, model_parameters, 3);

  bool store_output = false;
  if (parameters.size() > 4) {
    store_output = extract_bool_parameter(parameters, model_parameters, 4);
  }

  std::string scale_variable = "";
  if (parameters.size() > 5) {
    scale_variable = extract_string_parameter(parameters, model_parameters, 5);
  }

  arma::colvec scale = arma::colvec();
  if (parameters.size() > 6) {
    scale = extract_vector_parameter(parameters, model_parameters, 6);
  }

  bool parallel = false;
  if (parameters.size() > 7) {
    parallel = extract_bool_parameter(parameters, model_parameters, 7);
  }

  int grain_size = 100000;
  if (parameters.size() > 8) {
    grain_size = extract_int_parameter(parameters, model_parameters, 8);
  }

  return List::create(number_of_points, tolerance_variable, tolerance, p,
                      store_output, scale_variable, scale, parallel,
                      grain_size);
}

List get_abc_enki_parameter_info(const List &model_parameters,
                                 const List &current_sbi,
                                 const std::string &sbi_name) {

  // throw error if more than one variable
  // if (variable_names.size()!=1)
  //  Rcpp::stop("Only one variable allowed for " + sbi_name + ".");

  if (!current_sbi.containsElementNamed("parameters"))
    Rcpp::stop("Missing parameters for " + sbi_name +
               " (6-13 parameters required).");

  List parameters = current_sbi["parameters"];

  if (!((parameters.size() >= 6) && (parameters.size() <= 13)))
    Rcpp::stop("6-13 parameters required for " + sbi_name + ".");

  int number_of_points = extract_int_parameter(parameters, model_parameters, 0);

  int estimator_type_in =
      extract_int_parameter(parameters, model_parameters, 1);

  std::string tolerance_variable =
      extract_string_parameter(parameters, model_parameters, 2);

  double epsilon = extract_double_parameter(parameters, model_parameters, 3);

  std::string shifter_name =
      extract_string_parameter(parameters, model_parameters, 4);

  size_t number_of_targets_in =
      extract_int_parameter(parameters, model_parameters, 5);

  if (number_of_targets_in == 0) {
    Rcpp::stop("Error getting EnKI-ABC number of targets: must be at least 1.");
  }

  double significance_level;
  if (parameters.size() > 6) {
    significance_level =
        extract_double_parameter(parameters, model_parameters, 6);
  } else {
    significance_level = 1.0;
  }

  int enki_lag;
  if (parameters.size() > 7) {
    enki_lag = extract_int_parameter(parameters, model_parameters, 7);
  } else {
    enki_lag = 0;
  }

  std::string scale_variable_in;
  if (parameters.size() > 8) {
    scale_variable_in =
        extract_string_parameter(parameters, model_parameters, 8);
  } else {
    scale_variable_in = "";
  }

  arma::colvec scale_in;
  if (parameters.size() > 9) {
    scale_in = extract_vector_parameter(parameters, model_parameters, 9);
  } else {
    scale_in = arma::colvec();
  }

  bool enki_on_summary;

  if (parameters.size() > 10) {
    enki_on_summary = extract_bool_parameter(parameters, model_parameters, 10);
  } else {
    enki_on_summary = true;
  }

  bool parallel = false;
  if (parameters.size() > 11) {
    parallel = extract_bool_parameter(parameters, model_parameters, 11);
  }

  int grain_size = 100000;
  if (parameters.size() > 12) {
    grain_size = extract_int_parameter(parameters, model_parameters, 12);
  }

  return List::create(number_of_points, estimator_type_in, tolerance_variable,
                      epsilon, shifter_name, number_of_targets_in,
                      significance_level, enki_lag, scale_variable_in, scale_in,
                      enki_on_summary, parallel, grain_size);
}

List get_sl_parameter_info(const List &model_parameters,
                           const List &current_sbi,
                           const std::string &sbi_name) {

  if (!current_sbi.containsElementNamed("parameters"))
    Rcpp::stop("Missing parameters for " + sbi_name +
               " (1-5 parameters required).");

  List parameters = current_sbi["parameters"];

  if (!((parameters.size() == 1) || (parameters.size() == 2) ||
        (parameters.size() == 3) || (parameters.size() == 4) ||
        (parameters.size() == 5)))
    Rcpp::stop("1-5 parameters required for " + sbi_name + ".");

  int number_of_points = extract_int_parameter(parameters, model_parameters, 0);

  bool unbiased = false;
  if (parameters.size() > 1) {
    unbiased = extract_bool_parameter(parameters, model_parameters, 1);
  }

  bool store_output = false;
  if (parameters.size() > 2) {
    store_output = extract_bool_parameter(parameters, model_parameters, 2);
  }

  bool parallel = false;
  if (parameters.size() > 3) {
    parallel = extract_bool_parameter(parameters, model_parameters, 3);
  }

  int grain_size = 100000;
  if (parameters.size() > 4) {
    grain_size = extract_int_parameter(parameters, model_parameters, 4);
  }

  return List::create(number_of_points, unbiased, store_output, parallel,
                      grain_size);
}

List get_kf_info(const List &model_parameters, const List &current_likelihood,
                 const std::string &likelihood_name) {
  if (!current_likelihood.containsElementNamed("parameters"))
    Rcpp::stop("Missing parameters for " + likelihood_name +
               " (0-1 parameters required).");

  List parameters = current_likelihood["parameters"];

  if (!((parameters.size() == 0) || (parameters.size() == 1)))
    Rcpp::stop("0-1 parameters required for " + likelihood_name + ".");

  int iterations_to_store;
  if (parameters.size() > 0) {
    iterations_to_store =
        extract_int_parameter(parameters, model_parameters, 0);
  } else {
    iterations_to_store = 0;
  }

  return List::create(iterations_to_store);
}

List get_enkf_info(const List &model_parameters, const List &current_likelihood,
                   const std::string &likelihood_name) {
  if (!current_likelihood.containsElementNamed("parameters"))
    Rcpp::stop("Missing parameters for " + likelihood_name +
               " (1-4 parameters required).");

  List parameters = current_likelihood["parameters"];

  if (!((parameters.size() == 2) || (parameters.size() == 3) ||
        (parameters.size() == 4) || (parameters.size() == 5)))
    Rcpp::stop("2-5 parameters required for " + likelihood_name + ".");

  size_t number_of_ensemble_members =
      extract_int_parameter(parameters, model_parameters, 0);

  std::string shifter_name =
      extract_string_parameter(parameters, model_parameters, 1);

  int iterations_to_store;
  if (parameters.size() > 2) {
    iterations_to_store =
        extract_int_parameter(parameters, model_parameters, 2);
  } else {
    iterations_to_store = 0;
  }

  bool parallel;
  if (parameters.size() > 3) {
    parallel = extract_bool_parameter(parameters, model_parameters, 3);
  } else {
    parallel = false;
  }

  size_t grain_size;
  if (parameters.size() > 4) {
    grain_size = extract_int_parameter(parameters, model_parameters, 4);
  } else {
    grain_size = 100000;
  }

  return List::create(number_of_ensemble_members, shifter_name,
                      iterations_to_store, parallel, grain_size);
}

List get_pf_info(const List &model_parameters, const List &current_likelihood,
                 const std::string &likelihood_name) {
  if (!current_likelihood.containsElementNamed("parameters"))
    Rcpp::stop("Missing parameters for " + likelihood_name +
               " (1-4 parameters required).");

  List parameters = current_likelihood["parameters"];

  if (!((parameters.size() == 1) || (parameters.size() == 2) ||
        (parameters.size() == 3) || (parameters.size() == 4)))
    Rcpp::stop("1-4 parameters required for " + likelihood_name + ".");

  size_t number_of_particles =
      extract_int_parameter(parameters, model_parameters, 0);

  int iterations_to_store;
  if (parameters.size() > 1) {
    iterations_to_store =
        extract_int_parameter(parameters, model_parameters, 1);
  } else {
    // default: keep the final generation (its points and weights) for each
    // outer point (owner, 2026-11-15; was 0, i.e. nothing)
    iterations_to_store = 1;
  }

  bool parallel;
  if (parameters.size() > 2) {
    parallel = extract_bool_parameter(parameters, model_parameters, 2);
  } else {
    parallel = false;
  }

  size_t grain_size;
  if (parameters.size() > 3) {
    grain_size = extract_int_parameter(parameters, model_parameters, 3);
  } else {
    grain_size = 64; // as PF() (dev/plan.md 2.6b, G5)
  }

  return List::create(number_of_particles, iterations_to_store, parallel,
                      grain_size);
}

std::vector<size_t> get_enk_likelihood_indices(const List &model,
                                               const List &model_parameters) {
  std::vector<size_t> enk_likelihood_indices;
  if (model.containsElementNamed("method")) {
    List method_info = model["method"];

    for (size_t i = 0; i < method_info.size(); ++i) {
      List current_method = method_info[i];
      if (current_method.containsElementNamed("enk_likelihood_index")) {
        NumericVector indices;

        SEXP index_info_SEXP = current_method["enk_likelihood_index"];
        indices = load_numeric_vector(index_info_SEXP);

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

List get_ssm_info(const List &model, const List &model_parameters) {

  // throw error if more than one variable
  // if (index_names.size()!=1)
  //  Rcpp::stop("Only one index variable allowed for filter.");

  if (model.containsElementNamed("method")) {
    List method_info = model["method"];

    for (size_t i = 0; i < method_info.size(); ++i) {
      List current_method = method_info[i];
      if (current_method.containsElementNamed("ssm")) {
        List filtering_info = current_method["ssm"];
        if (filtering_info.containsElementNamed("parameters")) {
          if (Rf_isNewList(filtering_info["parameters"])) {
            List values = filtering_info["parameters"];

            std::string index_name_in =
                extract_string_parameter(values, model_parameters, 0);

            size_t first_index_in =
                extract_int_parameter(values, model_parameters, 1);

            size_t last_index_in =
                extract_int_parameter(values, model_parameters, 2);

            std::string time_name_in =
                extract_string_parameter(values, model_parameters, 3);

            double initial_time_in =
                extract_double_parameter(values, model_parameters, 4);

            std::string time_diff_name_in =
                extract_string_parameter(values, model_parameters, 5);

            double update_time_step_in =
                extract_double_parameter(values, model_parameters, 6);

            size_t predictions_per_update_in =
                extract_int_parameter(values, model_parameters, 7);

            return List::create(index_name_in, first_index_in, last_index_in,
                                time_name_in, initial_time_in,
                                time_diff_name_in, update_time_step_in,
                                predictions_per_update_in);
          } else {
            stop("Error in SSM specification.");
          }
        } else {
          stop("Error in SSM specification.");
        }
      }
    }
  } else {
    stop("Method not found in model.");
  }

  stop("SSM not found in model.");
}

bool is_linear_gaussian_data_model(const List &current_factor) {
  return current_factor.containsElementNamed("linear_gaussian_data_covariance") ||
         current_factor.containsElementNamed("linear_gaussian_data_model");
}

bool is_linear_gaussian_transition_model(const List &current_factor) {
  return current_factor.containsElementNamed(
             "linear_gaussian_transition_covariance") ||
         current_factor.containsElementNamed("linear_gaussian_transition_model");
}

// The one-line built-in: (first variable, [second variable,] matrix, covariance).
List get_one_line_linear_gaussian(const List &current_factor,
                                  const std::string &block_name,
                                  size_t number_of_variables,
                                  const List &model_parameters) {
  if (!Rf_isNewList(current_factor[block_name]))
    stop(block_name + ": invalid specification.");
  List info = current_factor[block_name];
  if (!(info.containsElementNamed("type") &&
        info.containsElementNamed("variables") &&
        info.containsElementNamed("parameters")))
    stop(block_name + ": invalid specification.");
  std::string type = info["type"];
  if (!iequals(type, block_name))
    stop(block_name + ": unknown type ilike::" + type + " (expected ilike::" +
         block_name + ").");
  std::string first_variable = info["variables"];
  List parameters = info["parameters"];
  size_t expected = number_of_variables + 1; // other variables, matrix, covariance
  if (parameters.size() != expected)
    stop("ilike::" + block_name + " requires " +
         std::to_string(expected + 1) + " arguments.");
  std::string second_variable = first_variable;
  size_t k = 0;
  if (number_of_variables == 2) {
    second_variable = Rcpp::as<std::string>(parameters[0]);
    k = 1;
  }
  arma::mat matrix = extract_matrix_parameter(parameters, model_parameters, k);
  arma::mat covariance =
      extract_matrix_parameter(parameters, model_parameters, k + 1);
  return List::create(first_variable, second_variable, matrix, covariance);
}

LinearGaussianSpec get_linear_gaussian_data_model(const List &current_factor,
                                                  const List &model_parameters) {
  LinearGaussianSpec spec;
  if (current_factor.containsElementNamed("linear_gaussian_data_model")) {
    List info = get_one_line_linear_gaussian(
        current_factor, "linear_gaussian_data_model", 2, model_parameters);
    spec.variable = Rcpp::as<std::string>(info[0]);
    spec.state_variable_string = Rcpp::as<std::string>(info[1]);
    spec.matrix = Rcpp::as<arma::mat>(info[2]);
    spec.covariance = Rcpp::as<arma::mat>(info[3]);
  } else {
    if (!current_factor.containsElementNamed("linear_gaussian_data_matrix"))
      stop("linear_gaussian_data must contain both matrix and covariance.");
    if (!current_factor.containsElementNamed("type"))
      stop("linear_gaussian_data must contain a type.");
    if (!current_factor.containsElementNamed("linear_gaussian_data_variable"))
      stop("linear_gaussian_data must contain linear_gaussian_data_variable.");
    if (!current_factor.containsElementNamed(
            "linear_gaussian_data_state_variable"))
      stop("linear_gaussian_data must contain "
           "linear_gaussian_data_state_variable.");
    SEXP matrix_SEXP = current_factor["linear_gaussian_data_matrix"];
    SEXP covariance_SEXP = current_factor["linear_gaussian_data_covariance"];
    spec.variable = load_string(current_factor["linear_gaussian_data_variable"]);
    spec.state_variable_string =
        load_string(current_factor["linear_gaussian_data_state_variable"]);
    size_t type = current_factor["type"];
    if (type == 1) {
      spec.matrix = load_matrix(matrix_SEXP);
      spec.covariance = load_matrix(covariance_SEXP);
    } else if (type == 2) {
      spec.is_function = true;
      spec.matrix_function = load_matrix_function(matrix_SEXP);
      spec.covariance_function = load_matrix_function(covariance_SEXP);
    } else {
      stop("linear_gaussian_data specificiation is invalid.");
    }
  }
  spec.state_variables = split(spec.state_variable_string, ';');
  return spec;
}

LinearGaussianSpec
get_linear_gaussian_transition_model(const List &current_factor,
                                     const List &model_parameters,
                                     bool variable_required) {
  LinearGaussianSpec spec;
  if (current_factor.containsElementNamed("linear_gaussian_transition_model")) {
    List info = get_one_line_linear_gaussian(
        current_factor, "linear_gaussian_transition_model", 1,
        model_parameters);
    spec.variable = Rcpp::as<std::string>(info[0]);
    spec.matrix = Rcpp::as<arma::mat>(info[2]);
    spec.covariance = Rcpp::as<arma::mat>(info[3]);
  } else {
    if (!current_factor.containsElementNamed("linear_gaussian_transition_matrix"))
      stop("linear_gaussian must contain both matrix and covariance.");
    if (!current_factor.containsElementNamed("type"))
      stop("linear_gaussian must contain a type.");
    if (current_factor.containsElementNamed(
            "linear_gaussian_transition_variable"))
      spec.variable =
          load_string(current_factor["linear_gaussian_transition_variable"]);
    else if (variable_required)
      stop("linear_gaussian must contain linear_gaussian_transition_variable.");
    SEXP matrix_SEXP = current_factor["linear_gaussian_transition_matrix"];
    SEXP covariance_SEXP =
        current_factor["linear_gaussian_transition_covariance"];
    size_t type = current_factor["type"];
    if (type == 1) {
      spec.matrix = load_matrix(matrix_SEXP);
      spec.covariance = load_matrix(covariance_SEXP);
    } else if (type == 2) {
      spec.is_function = true;
      spec.matrix_function = load_matrix_function(matrix_SEXP);
      spec.covariance_function = load_matrix_function(covariance_SEXP);
    } else {
      stop("linear_gaussian specificiation is invalid.");
    }
  }
  spec.state_variable_string = spec.variable;
  spec.state_variables = std::vector<std::string>(1, spec.variable);
  return spec;
}
