#include <RcppArmadillo.h>
using namespace Rcpp;

#include <string>
#include <sstream>
#include <vector>
#include <chrono>

#include "custom_no_params_proposal_kernel.h"
#include "custom_proposal_kernel.h"
#include "distributions.h"
#include "parameters.h"
#include "exact_likelihood_estimator.h"
#include "importance_sampler.h"
#include "smc_mcmc_move.h"
#include "custom_distribution_factor.h"
#include "gaussian_distribution_factor.h"
#include "loggaussian_distribution_factor.h"
#include "gamma_distribution_factor.h"
#include "custom_likelihood_factor.h"
#include "gaussian_independent_proposal_kernel.h"
#include "loggaussian_independent_proposal_kernel.h"
#include "gamma_independent_proposal_kernel.h"
#include "custom_distribution_proposal_kernel.h"
#include "custom_independent_proposal_kernel.h"
#include "custom_guided_distribution_proposal_kernel.h"
#include "custom_guided_independent_proposal_kernel.h"
#include "composite_independent_proposal_kernel.h"
#include "smc_output.h"
#include "metropolis_mcmc.h"
#include "metropolis_hastings_mcmc.h"
#include "unadjusted_mcmc.h"
#include "stochastic_scan_mcmc.h"
#include "deterministic_scan_mcmc.h"
#include "gaussian_random_walk_proposal_kernel.h"
#include "uniform_random_walk_proposal_kernel.h"
#include "langevin_proposal_kernel.h"
#include "barker_dynamics_proposal_kernel.h"
#include "hmc_proposal_kernel.h"
#include "mirror_proposal_kernel.h"
#include "custom_guided_no_params_proposal_kernel.h"
#include "custom_no_params_proposal_kernel.h"
#include "custom_proposal_kernel.h"
#include "custom_guided_proposal_kernel.h"
#include "custom_guided_no_params_symmetric_proposal_kernel.h"
#include "custom_no_params_symmetric_proposal_kernel.h"
#include "custom_symmetric_proposal_kernel.h"
#include "custom_guided_symmetric_proposal_kernel.h"
#include "iterations_mcmc_termination.h"
#include "se_mcmc_termination.h"
#include "ess_smc_criterion.h"
#include "cess_smc_criterion.h"
#include "positive_smc_criterion.h"
#include "utils.h"
#include "annealed_likelihood_estimator.h"
#include "stable_smc_termination.h"
#include "always_smc_termination.h"
#include "direct_gradient_estimator.h"
#include "vector_single_index.h"
#include "likelihood_maker.h"
#include "transform.h"
#include "measurement_covariance_estimator.h"
#include "ensemble_kalman_inversion.h"
#include "ensemble_kalman_output.h"
#include "linear_gaussian_noise_proposal_kernel.h"
#include "linear_gaussian_noise_function_proposal_kernel.h"
#include "nonlinear_gaussian_noise_proposal_kernel.h"
#include "nonlinear_gaussian_noise_function_proposal_kernel.h"
#include "generic_measurement_covariance_estimator.h"
#include "direct_linear_gaussian_measurement_covariance_estimator.h"
#include "direct_nonlinear_gaussian_measurement_covariance_estimator.h"
#include "mixed_generic_direct_gaussian_measurement_covariance_estimator.h"
#include "stochastic_ensemble_shifter.h"
#include "kalman_filter.h"
#include "kalman_filter_output.h"
#include "exact_kalman_predictor.h"
#include "exact_kalman_updater.h"

#include "linear_gaussian_state_space_model.h"

// from https://stackoverflow.com/questions/2165921/converting-from-a-stdstring-to-bool
bool stob(const std::string &s)
{
  auto result = false;    // failure to assert is false
  
  std::istringstream is(s);
  // first try simple integer conversion
  is >> result;
  
  if (is.fail())
  {
    // simple integer failed; try boolean
    is.clear();
    is >> std::boolalpha >> result;
  }
  
  if (is.fail())
  {
    stop(s + " is not convertable to bool");
  }
  
  return result;
}

std::vector<std::string> split(const std::string &str, const char delimiter)
{
  std::vector<std::string> tokens;
  std::stringstream ss(str);
  std::string token;
  while (std::getline(ss, token, delimiter)) {
    tokens.push_back(token);
  }
  return tokens;
}

bool isDoubleInList(const Rcpp::List &list, int index_from_one)
{
  if (index_from_one<=list.size())
  {
    SEXP var = list[index_from_one-1];
    return Rf_isReal(var) && Rf_length(var) == 1;
  }
  else
  {
    return false;
  }
}

bool isStringInList(const Rcpp::List &list, int index_from_one)
{
  if (index_from_one<=list.size())
  {
    SEXP var = list[index_from_one-1];
    return Rf_isString(var) && Rf_length(var) == 1;
  }
  else
  {
    return false;
  }
}

bool isVectorInList(const Rcpp::List &list, int index_from_one)
{
  if (index_from_one<=list.size())
  {
    SEXP var = list[index_from_one-1];
    return Rf_isVector(var);
  }
  else
  {
    return false;
  }
}

bool isMatrixInList(const Rcpp::List &list, int index_from_one)
{
  if (index_from_one<=list.size())
  {
    SEXP var = list[index_from_one-1];
    return Rf_isMatrix(var);
  }
  else
  {
    return false;
  }
}

bool isNumericVectorInList(const Rcpp::List &list, int index_from_one)
{
  if (index_from_one<=list.size())
  {
    SEXP var = list[index_from_one-1];
    return Rf_isVector(var) && Rf_isNumeric(var);
  }
  else
  {
    return false;
  }
}

bool isNumericMatrixInList(const Rcpp::List &list, int index_from_one)
{
  if (index_from_one<=list.size())
  {
    SEXP var = list[index_from_one-1];
    return Rf_isMatrix(var) && Rf_isNumeric(var);
  }
  else
  {
    return false;
  }
}

bool isStringVectorInList(const Rcpp::List &list, const std::string &name)
{
  SEXP var = list[name];
  return Rf_isVector(var) && Rf_isString(var);
}

bool isDoubleInList(const Rcpp::List &list, const std::string &name)
{
  SEXP var = list[name];
  return Rf_isReal(var) && Rf_length(var)==1;
}

bool isNumericVectorInList(const Rcpp::List &list, const std::string &name)
{
  SEXP var = list[name];
  return Rf_isVector(var) && Rf_isNumeric(var);
}

bool isNumericMatrixInList(const Rcpp::List &list, const std::string &name)
{
  SEXP var = list[name];
  return Rf_isMatrix(var) && Rf_isNumeric(var);
}

/*
Data example_data()
{
  RandomNumberGenerator rng;
  rng.seed(1);
  size_t n = 100;
  //arma::colvec sampled = rnorm(rng,n);
  arma::colvec sampled(n);
  for (size_t i = 0; i<n; ++i)
    sampled[i] = rnorm(rng,0.0,1.0);
  Data data;
  data["y"] = sampled;
  return data;
}

XPtr<DataPtr> example_store_data();

// [[Rcpp::export]]
XPtr<DataPtr> example_store_data()
{
  return(XPtr<DataPtr>(new DataPtr(&example_data)));
}
*/

Data get_data(const List &model)
{
  if (model.containsElementNamed("data"))
  {
    List data_list = model["data"];
    List data_list0 = data_list[0];
    SEXP data_SEXP = data_list0["data"];
    
    Data output = load_data(data_SEXP);
    return output;
    
    //return load_data(data_SEXP);
  }
  else
  {
    Rcpp::stop("do_importance_sampler: data not found in model specification.");
  }
}

bool extract_bool_parameter(const List &parameters_from_file,
                            const List &model_parameters,
                            size_t index)
{
  bool result;
  std::string parameter_string = parameters_from_file[index];
  
  if (parameter_string.at(0)=='p')
  {
    size_t parameter_index;
    try
    {
      parameter_index = std::stoi(parameter_string.substr(1));
    } catch (...)
    {
      Rcpp::stop("Parameter index in model file not an integer.");
    }
    
    if (isDoubleInList(model_parameters,parameter_index))
    {
      result = model_parameters[parameter_index-1];
    }
    else
    {
      Rcpp::stop("extract_bool_parameter: Parameter index does not correspond to a real number in the parameters file.");
    }
    
  }
  else
  {
    try
    {
      result = stob(parameter_string);
    }
    catch (...)
    {
      Rcpp::stop("Parameter in model file is not a real number.");
    }
  }
  
  return result;
}

int extract_int_parameter(const List &parameters_from_file,
                          const List &model_parameters,
                          size_t index)
{
  int result;
  std::string parameter_string = parameters_from_file[index];
  
  if (parameter_string.at(0)=='p')
  {
    size_t parameter_index;
    try
    {
      parameter_index = std::stoi(parameter_string.substr(1));
    } catch (...)
    {
      Rcpp::stop("Parameter index in model file not an integer.");
    }
    
    if (isDoubleInList(model_parameters,parameter_index))
    {
      result = model_parameters[parameter_index-1];
    }
    else
    {
      Rcpp::stop("extract_int_parameter: Parameter index does not correspond to a real number in the parameters file. Did you supply a model_parameter_list argument?");
    }
    
  }
  else
  {
    try
    {
      result = std::stoi(parameter_string);
    }
    catch (...)
    {
      Rcpp::stop("Parameter in model file is not a real number.");
    }
  }
  
  return result;
}

double extract_double_parameter(const List &parameters_from_file,
                                const List &model_parameters,
                                size_t index)
{
  double result;
  std::string parameter_string = parameters_from_file[index];
  
  if (parameter_string.at(0)=='p')
  {
    size_t parameter_index;
    try
    {
      parameter_index = std::stoi(parameter_string.substr(1));
    } catch (...)
    {
      Rcpp::stop("Parameter index in model file not an integer.");
    }
    
    if (isDoubleInList(model_parameters,parameter_index))
    {
      result = model_parameters[parameter_index-1];
    }
    else
    {
      Rcpp::stop("extract_double_parameter: Parameter index does not correspond to a real number in the parameters file.");
    }
    
  }
  else
  {
    try
    {
      result = std::stod(parameter_string);
    }
    catch (...)
    {
      Rcpp::stop("Parameter in model file is not a real number.");
    }
  }
  
  return result;
}

std::string extract_string_parameter(const List &parameters_from_file,
                                     const List &model_parameters,
                                     size_t index)
{
  std::string parameter_string = parameters_from_file[index];
  return parameter_string;
}

arma::colvec extract_vector_parameter(const List &parameters_from_file,
                                      const List &model_parameters,
                                      size_t index)
{
  arma::colvec result;
  std::string parameter_string = parameters_from_file[index];
  
  if (parameter_string.at(0)=='p')
  {
    size_t parameter_index;
    try
    {
      parameter_index = std::stoi(parameter_string.substr(1));
    } catch (...)
    {
      Rcpp::stop("Parameter index in model file not an integer.");
    }
    
    if (isVectorInList(model_parameters,parameter_index))
    {
      result = as<NumericVector>(model_parameters[parameter_index-1]);
    }
    else
    {
      Rcpp::stop("extract_vector_parameter: Parameter index does not correspond to a real number in the parameters file.");
    }
    
  }
  else
  {
    try
    {
      result = std::stod(parameter_string);
    }
    catch (...)
    {
      Rcpp::stop("Parameter in model file is not a real number.");
    }
  }
  
  return result;
}

arma::mat extract_matrix_parameter(const List &parameters_from_file,
                                   const List &model_parameters,
                                   size_t index)
{
  arma::mat result;
  std::string parameter_string = parameters_from_file[index];
  
  if (parameter_string.at(0)=='p')
  {
    size_t parameter_index;
    try
    {
      parameter_index = std::stoi(parameter_string.substr(1));
    } catch (...)
    {
      Rcpp::stop("Parameter index in model file not an integer.");
    }
    
    if (isMatrixInList(model_parameters,parameter_index))
    {
      result = as<arma::mat>(model_parameters[parameter_index-1]);
    }
    else
    {
      Rcpp::stop("extract_matrix_parameter: Parameter index does not correspond to a real number in the parameters file.");
    }
    
  }
  else
  {
    try
    {
      result = std::stod(parameter_string);
    }
    catch (...)
    {
      Rcpp::stop("Parameter in model file is not a real number.");
    }
  }
  
  return result;
}

List get_single_variable_one_double_parameter_info(const List &model_parameters,
                                                   const List &current_distribution,
                                                   const std::string &distribution_name)
{
  std::string augmented_variable_names = current_distribution["variables"];
  
  // split string in ;
  std::vector<std::string> variable_names = split(augmented_variable_names,';');
  
  // throw error if more than one variable
  if (variable_names.size()!=1)
    Rcpp::stop("Only one variable allowed for " + distribution_name + ".");
  
  if (!current_distribution.containsElementNamed("parameters"))
    Rcpp::stop("Missing parameters for " + distribution_name + " (one parameters required).");
  
  List parameters = current_distribution["parameters"];
  
  if (parameters.size()!=1)
    Rcpp::stop("One parameter required for " + distribution_name + ".");
  
  double first_param = extract_double_parameter(parameters,
                                                model_parameters,
                                                0);
  
  return List::create(variable_names[0],first_param);
}

List get_single_variable_two_double_parameter_info(const List &model_parameters,
                                                   const List &current_distribution,
                                                   const std::string &distribution_name)
{
  std::string augmented_variable_names = current_distribution["variables"];
  
  // split string in ;
  std::vector<std::string> variable_names = split(augmented_variable_names,';');
  
  // throw error if more than one variable
  if (variable_names.size()!=1)
    Rcpp::stop("Only one variable allowed for " + distribution_name + ".");
  
  if (!current_distribution.containsElementNamed("parameters"))
    Rcpp::stop("Missing parameters for " + distribution_name + " (two parameters required).");
  
  List parameters = current_distribution["parameters"];
  
  if (parameters.size()!=2)
    Rcpp::stop("Two parameters required for " + distribution_name + ".");
  
  double first_param = extract_double_parameter(parameters,
                                                model_parameters,
                                                0);
  
  double second_param = extract_double_parameter(parameters,
                                                 model_parameters,
                                                 1);
  
  return List::create(variable_names[0],first_param,second_param);
}

List get_single_variable_vector_parameter_info(const List &model_parameters,
                                               const List &current_distribution,
                                               const std::string &distribution_name)
{
  std::string augmented_variable_names = current_distribution["variables"];
  
  // split string in ;
  std::vector<std::string> variable_names = split(augmented_variable_names,';');
  
  // throw error if more than one variable
  if (variable_names.size()!=1)
    Rcpp::stop("Only one variable allowed for " + distribution_name + ".");
  
  if (!current_distribution.containsElementNamed("parameters"))
    Rcpp::stop("Missing parameters for " + distribution_name + " (one parameter required).");
  
  List parameters = current_distribution["parameters"];
  
  if (parameters.size()!=1)
    Rcpp::stop("One parameter required for " + distribution_name + ".");
  
  arma::colvec first_param = extract_vector_parameter(parameters,
                                                      model_parameters,
                                                      0);
  
  return List::create(variable_names[0],first_param);
}

List get_single_variable_matrix_parameter_info(const List &model_parameters,
                                               const List &current_distribution,
                                               const std::string &distribution_name)
{
  std::string augmented_variable_names = current_distribution["variables"];
  
  // split string in ;
  std::vector<std::string> variable_names = split(augmented_variable_names,';');
  
  // throw error if more than one variable
  if (variable_names.size()!=1)
    Rcpp::stop("Only one variable allowed for " + distribution_name + ".");
  
  if (!current_distribution.containsElementNamed("parameters"))
    Rcpp::stop("Missing parameters for " + distribution_name + " (one parameter required).");
  
  List parameters = current_distribution["parameters"];
  
  if (parameters.size()!=1)
    Rcpp::stop("One parameter required for " + distribution_name + ".");
  
  arma::mat first_param = extract_matrix_parameter(parameters,
                                                   model_parameters,
                                                   0);
  
  return List::create(variable_names[0],first_param);
}

List get_single_variable_vector_and_matrix_parameter_info(const List &model_parameters,
                                                          const List &current_distribution,
                                                          const std::string &distribution_name)
{
  std::string augmented_variable_names = current_distribution["variables"];
  
  // split string in ;
  std::vector<std::string> variable_names = split(augmented_variable_names,';');
  
  // throw error if more than one variable
  if (variable_names.size()!=1)
    Rcpp::stop("Only one variable allowed for " + distribution_name + ".");
  
  if (!current_distribution.containsElementNamed("parameters"))
    Rcpp::stop("Missing parameters for " + distribution_name + " (two parameters required).");
  
  List parameters = current_distribution["parameters"];
  
  if (parameters.size()!=2)
    Rcpp::stop("Two parameters required for " + distribution_name + ".");
  
  arma::colvec first_param = extract_vector_parameter(parameters,
                                                      model_parameters,
                                                      0);
  
  arma::mat second_param = extract_matrix_parameter(parameters,
                                                    model_parameters,
                                                    1);
  
  return List::create(variable_names[0],first_param,second_param);
}

List get_abc_euclidean_uniform_parameter_info(const List &model_parameters,
                                              const List &current_sbi,
                                              const std::string &sbi_name)
{
  std::string augmented_variable_names = current_sbi["variables"];
  
  // split string in ;
  std::vector<std::string> variable_names = split(augmented_variable_names,';');
  
  if (!current_sbi.containsElementNamed("parameters"))
    Rcpp::stop("Missing parameters for " + sbi_name + " (3-5 parameters required).");
  
  List parameters = current_sbi["parameters"];
  
  if (! ( (parameters.size()==3) || (parameters.size()==4) || (parameters.size()==5) ) )
    Rcpp::stop("3-5 parameters required for " + sbi_name + ".");
  
  int number_of_points = extract_int_parameter(parameters,
                                               model_parameters,
                                               0);
  
  std::string tolerance_variable = extract_string_parameter(parameters,
                                                            model_parameters,
                                                            1);
  
  double tolerance = extract_double_parameter(parameters,
                                              model_parameters,
                                              2);
  
  bool parallel = false;
  if (parameters.size()>3)
  {
    parallel = extract_bool_parameter(parameters,
                                      model_parameters,
                                      3);
  }
  
  int grain_size = 100000;
  if (parameters.size()>4)
  {
    grain_size = extract_int_parameter(parameters,
                                       model_parameters,
                                       4);
  }
  
  return List::create(variable_names,number_of_points,tolerance_variable,tolerance,parallel,grain_size);
}

List get_abc_lp_uniform_parameter_info(const List &model_parameters,
                                       const List &current_sbi,
                                       const std::string &sbi_name)
{
  std::string augmented_variable_names = current_sbi["variables"];
  
  // split string in ;
  std::vector<std::string> variable_names = split(augmented_variable_names,';');
  
  if (!current_sbi.containsElementNamed("parameters"))
    Rcpp::stop("Missing parameters for " + sbi_name + " (4-6 parameters required).");
  
  List parameters = current_sbi["parameters"];
  
  if (! ( (parameters.size()==4) || (parameters.size()==5) || (parameters.size()==6) ) )
    Rcpp::stop("4-6 parameters required for " + sbi_name + ".");
  
  int number_of_points = extract_int_parameter(parameters,
                                               model_parameters,
                                               0);
  
  std::string tolerance_variable = extract_string_parameter(parameters,
                                                            model_parameters,
                                                            1);
  
  double tolerance = extract_double_parameter(parameters,
                                              model_parameters,
                                              2);
  
  double p = extract_double_parameter(parameters,
                                      model_parameters,
                                      3);
  
  bool parallel = false;
  if (parameters.size()>4)
  {
    parallel = extract_bool_parameter(parameters,
                                      model_parameters,
                                      4);
  }
  
  int grain_size = 100000;
  if (parameters.size()>5)
  {
    grain_size = extract_int_parameter(parameters,
                                       model_parameters,
                                       5);
  }
  
  return List::create(variable_names,number_of_points,tolerance_variable,tolerance,p,parallel,grain_size);
}

List get_abc_enki_parameter_info(const List &model_parameters,
                                 const List &current_sbi,
                                 const std::string &sbi_name)
{
  std::string augmented_variable_names = current_sbi["variables"];
  
  // split string in ;
  std::vector<std::string> variable_names = split(augmented_variable_names,';');
  
  // throw error if more than one variable
  //if (variable_names.size()!=1)
  //  Rcpp::stop("Only one variable allowed for " + sbi_name + ".");
  
  if (!current_sbi.containsElementNamed("parameters"))
    Rcpp::stop("Missing parameters for " + sbi_name + " (8-10 parameters required).");
  
  List parameters = current_sbi["parameters"];
  
  if (! ( (parameters.size()==8) || (parameters.size()==9) || (parameters.size()==10) ) )
    Rcpp::stop("8-10 parameters required for " + sbi_name + ".");
  
  int number_of_points = extract_int_parameter(parameters,
                                               model_parameters,
                                               0);
  Rcout << number_of_points << std::endl;
  
  std::string tolerance_variable = extract_string_parameter(parameters,
                                                            model_parameters,
                                                            1);
  Rcout << number_of_points << std::endl;
  arma::colvec schedule_colvec = extract_vector_parameter(parameters,
                                                          model_parameters,
                                                          2);
  Rcout << number_of_points << std::endl;
  
  std::vector<double> schedule;
  if (schedule_colvec.n_elem==0)
  {
    stop("Error getting EnKI-ABC schedule: need at least one tolerance value.");
  }
  else if (schedule_colvec.n_elem==1)
  {
    schedule.push_back(arma::datum::inf);
    schedule.push_back(schedule_colvec[0]);
  }
  else
  {
    for (size_t k=0; k<schedule_colvec.n_elem; ++k)
    {
      schedule.push_back(schedule_colvec[k]);
    }
  }
  
  int enki_lag = extract_int_parameter(parameters,
                                       model_parameters,
                                       3);
  Rcout << number_of_points << std::endl;
  std::string shifter_name = extract_string_parameter(parameters,
                                                      model_parameters,
                                                      4);
  Rcout << number_of_points << std::endl;
  
  double proportion = extract_double_parameter(parameters,
                                               model_parameters,
                                               5);
  Rcout << number_of_points << std::endl;
  
  double enki_annealing_desired_cess = proportion*double(number_of_points);
  
  int enki_number_of_bisections = extract_int_parameter(parameters,
                                                        model_parameters,
                                                        6);
  Rcout << number_of_points << std::endl;
  
  bool enki_on_summary = extract_bool_parameter(parameters,
                                                model_parameters,
                                                7);
  Rcout << number_of_points << std::endl;
  
  double significance_level = extract_double_parameter(parameters,
                                                       model_parameters,
                                                       8);
  Rcout << number_of_points << std::endl;
    
  bool parallel = false;
  if (parameters.size()>9)
  {
    parallel = extract_bool_parameter(parameters,
                                      model_parameters,
                                      9);
  }
  
  int grain_size = 100000;
  if (parameters.size()>10)
  {
    grain_size = extract_int_parameter(parameters,
                                       model_parameters,
                                       10);
  }
  
  Rcout << number_of_points << std::endl;
  return List::create(variable_names,
                      number_of_points,
                      tolerance_variable,
                      schedule,
                      enki_lag,
                      shifter_name,
                      enki_annealing_desired_cess,
                      enki_number_of_bisections,
                      enki_on_summary,
                      significance_level,
                      parallel,
                      grain_size);
}


List get_filtering_info(const List &model_parameters,
                        const List &filtering_info)
{
  //std::string index_name_in = filtering_info["variables"];
  
  // split string in ;
  //std::vector<std::string> index_names = split(index_name_in,';');
  
  // throw error if more than one variable
  //if (index_names.size()!=1)
  //  Rcpp::stop("Only one index variable allowed for filter.");
  
  if (filtering_info.containsElementNamed("values"))
  {
    if (Rf_isNewList(filtering_info["values"]))
    {
      List values = filtering_info["values"];
      
      std::string index_name_in = extract_string_parameter(values,
                                                           model_parameters,
                                                           0);
      //output[0] = index_name_in;
      //Rcout << index_name_in << std::endl;
      size_t first_index_in = extract_int_parameter(values,
                                                    model_parameters,
                                                    1);
      //Rcout << first_index_in << std::endl;
      //output[1] = first_index_in;
      size_t last_index_in = extract_int_parameter(values,
                                                   model_parameters,
                                                   2);
      //Rcout << last_index_in << std::endl;
      //output[2] = last_index_in;
      std::string time_name_in = extract_string_parameter(values,
                                                          model_parameters,
                                                          3);
      //Rcout << time_name_in << std::endl;
      //output[3] = time_name_in;
      double initial_time_in = extract_double_parameter(values,
                                                        model_parameters,
                                                        4);
      //Rcout << initial_time_in << std::endl;
      //output[4] = initial_time_in;
      std::string time_diff_name_in = extract_string_parameter(values,
                                                               model_parameters,
                                                               5);
      //Rcout << time_diff_name_in << std::endl;
      //output[5] = time_diff_name_in;
      double update_time_step_in = extract_double_parameter(values,
                                                            model_parameters,
                                                            6);
      //Rcout << update_time_step_in << std::endl;
      //output[6] = update_time_step_in;
      size_t predictions_per_update_in = extract_int_parameter(values,
                                                               model_parameters,
                                                               7);
      //Rcout << predictions_per_update_in << std::endl;
      //output[7] = predictions_per_update_in;
      std::string state_name_in = extract_string_parameter(values,
                                                           model_parameters,
                                                           8);
      //Rcout << state_name_in << std::endl;
      //output[8] = state_name_in;
      std::string measurement_name_in = extract_string_parameter(values,
                                                                 model_parameters,
                                                                 9);
      //Rcout << measurement_name_in << std::endl;
      //output[9] = measurement_name_in;
      
      return List::create(index_name_in,
                          first_index_in,
                          last_index_in,
                          time_name_in,
                          initial_time_in,
                          time_diff_name_in,
                          update_time_step_in,
                          predictions_per_update_in,
                          state_name_in,
                          measurement_name_in);
    }
    else
    {
      stop("Error in filter method.");
    }
  }
  else
  {
    stop("Error in filter method.");
  }
}

std::vector<LikelihoodEstimator*> get_likelihood_estimators(RandomNumberGenerator* rng_in,
                                                            size_t* seed_in,
                                                            Data* data_in,
                                                            const List &model,
                                                            const List &model_parameters,
                                                            bool include_priors,
                                                            const std::vector<std::string> &sequencer_types,
                                                            const std::vector<std::string> &sequencer_variables,
                                                            const std::vector<std::vector<double>> &sequencer_schedules,
                                                            IndependentProposalKernel* proposal_in,
                                                            Index* &without_cancelled_index,
                                                            Index* &full_index,
                                                            bool &any_annealing,
                                                            const std::vector<int> &factors_affected_by_smc_sequence,
                                                            std::vector<Data> &data_created_in_get_likelihood_estimators)
{
  data_created_in_get_likelihood_estimators.clear();
  
  std::vector<size_t> without_cancelled_index_vector;
  std::vector<size_t> full_index_vector;
  
  any_annealing = false;
  std::string annealing_variable;
  for (size_t i=0; i<sequencer_types.size(); ++i)
  {
    //Rcout << sequencer_types[i] << std::endl;
    //Rcout << sequencer_schedules[i][0] << std::endl;
    //Rcout << sequencer_schedules[i][1] << std::endl;
    //if ( ( (sequencer_types[i]=="annealing") || (sequencer_types[i]=="tempering") ) && (!( (sequencer_schedules[i].size()==2) && (sequencer_schedules[i][0]==0.0) && (sequencer_schedules[i][1]==1.0)) ) )
    if ( (sequencer_types[i]=="annealing") || (sequencer_types[i]=="tempering") )
    {
      if (any_annealing==false)
      {
        any_annealing = true;
        annealing_variable = sequencer_variables[i];
        break;
      }
      else
      {
        stop("Two annealing variables specified: not currently supported.");
      }
    }
  }
  
  PowerFunctionPtr power = annealing_power;
  PowerFunctionPtr second_power = annealing_one_minus_power;
  
  std::vector<LikelihoodEstimator*> likelihood_estimators;
  //std::vector<DistributionFactor*> numerator_distribution_factors;
  //std::vector<LikelihoodFactor*> numerator_likelihood_factors;
  
  //std::vector<DistributionFactor*> prior_factors;
  //std::vector<LikelihoodFactor*> exact_likelihood_factors;
  
  //std::vector<IndependentProposalKernel*> numerator_proposal_factors;

  if ( model.containsElementNamed("factor") )
  {
    
    List factors = model["factor"];

    for (size_t i=0; i<factors.size(); ++i)
    {
      
      bool factor_is_fixed_in_smc = true;
      for (size_t j=0; j<factors_affected_by_smc_sequence.size(); ++j)
      {
        if (i+1==factors_affected_by_smc_sequence[j])
          factor_is_fixed_in_smc = false;
      }

      if (Rf_isNewList(factors[i]))
      {
        
        List current_factor = factors[i];
        if ( current_factor.containsElementNamed("prior") )
        {
          DistributionFactor* new_factor;
          
          if (Rf_isNewList(current_factor["prior"]))
          {
            
            List current_prior = current_factor["prior"];
            if ( current_prior.containsElementNamed("type") && current_prior.containsElementNamed("variables") )
            {
              std::string type = current_prior["type"];
              
              
              if (type=="norm")
              {
                List info = get_single_variable_two_double_parameter_info(model_parameters,
                                                                          current_prior,
                                                                          type);
                
                std::string variable = info[0];
                double mean = info[1];
                double sd = info[2];
                
                new_factor = new GaussianDistributionFactor(variable,
                                                            mean,
                                                            sd);
                
              }
              else if (type=="mvnorm")
              {
                List info = get_single_variable_vector_and_matrix_parameter_info(model_parameters,
                                                                                 current_prior,
                                                                                 type);
                
                std::string variable = info[0];
                arma::colvec mean = info[1];
                arma::mat cov = info[2];
                
                new_factor = new GaussianDistributionFactor(variable,
                                                            mean,
                                                            cov);
                
              }
              else if (type=="lnorm")
              {
                List info = get_single_variable_two_double_parameter_info(model_parameters,
                                                                          current_prior,
                                                                          type);
                
                std::string variable = info[0];
                double mean = info[1];
                double sd = info[2];
                
                new_factor = new LogGaussianDistributionFactor(variable,
                                                               mean,
                                                               sd);
                
              }
              else if (type=="mvlnorm")
              {
                List info = get_single_variable_vector_and_matrix_parameter_info(model_parameters,
                                                                                 current_prior,
                                                                                 type);
                
                std::string variable = info[0];
                arma::colvec mean = info[1];
                arma::mat cov = info[2];
                
                new_factor = new LogGaussianDistributionFactor(variable,
                                                               mean,
                                                               cov);
                
              }
              else if (type=="gamma")
              {
                List info = get_single_variable_two_double_parameter_info(model_parameters,
                                                                          current_prior,
                                                                          type);
                
                std::string variable = info[0];
                double shape = info[1];
                double rate = info[2];
                
                new_factor = new GammaDistributionFactor(variable,
                                                         shape,
                                                         rate);
                
              }
              else
              {
                Rcout << "Prior type " << type;
                stop("Prior type unknown.");
              }
              
              LikelihoodEstimator* new_likelihood_estimator = new ExactLikelihoodEstimator(rng_in,
                                                                                           seed_in,
                                                                                           data_in,
                                                                                           new_factor,
                                                                                           true);
              
              
              full_index_vector.push_back(likelihood_estimators.size());
              if (include_priors)
              {
                without_cancelled_index_vector.push_back(likelihood_estimators.size());
                
                if (any_annealing)
                {
                  likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(rng_in,
                                                                                  seed_in,
                                                                                  data_in,
                                                                                  new_likelihood_estimator,
                                                                                  power,
                                                                                  annealing_variable,
                                                                                  false));
                }
                else
                {
                  likelihood_estimators.push_back(new_likelihood_estimator);
                }
              }
              else
              {
                likelihood_estimators.push_back(new_likelihood_estimator);
              }
              
            }
            else
            {
              stop("Missing information for prior in model file.");
            }
            
          }
          else
          {
            stop("Error in factors section of model file.");
          }
          
        }
        else if (current_factor.containsElementNamed("evaluate_log_prior"))
        {
          SEXP evaluate_log_prior_SEXP = current_factor["evaluate_log_prior"];
          DistributionFactor* new_factor = new CustomDistributionFactor(load_evaluate_log_distribution(evaluate_log_prior_SEXP));
          LikelihoodEstimator* new_likelihood_estimator = new ExactLikelihoodEstimator(rng_in,
                                                                                       seed_in,
                                                                                       data_in,
                                                                                       new_factor,
                                                                                       factor_is_fixed_in_smc);
          full_index_vector.push_back(likelihood_estimators.size());
                             
          //Rcout << "loading evaluate_log_prior" << std::endl;
          if (include_priors)
          {
            //Rcout << "include priors" << std::endl;
            
            without_cancelled_index_vector.push_back(likelihood_estimators.size());
            
            if (any_annealing)
            {
              //Rcout << "any annealing" << std::endl;
              likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(rng_in,
                                                                              seed_in,
                                                                              data_in,
                                                                              new_likelihood_estimator,
                                                                              power,
                                                                              annealing_variable,
                                                                              false));
            }
            else
            {
              //Rcout << "no annealing" << std::endl;
              likelihood_estimators.push_back(new_likelihood_estimator);
            }
          }
          else
          {
            //Rcout << "not including priors" << std::endl;
            likelihood_estimators.push_back(new_likelihood_estimator);
          }
          
        }
        else if ( current_factor.containsElementNamed("simulate_prior") )
        {
          
        }
        else if ( current_factor.containsElementNamed("evaluate_log_likelihood") )
        {
          
          SEXP evaluate_log_likelihood_SEXP = current_factor["evaluate_log_likelihood"];
          LikelihoodFactor* new_factor = new CustomLikelihoodFactor(load_evaluate_log_likelihood(evaluate_log_likelihood_SEXP),
                                                                    data_in);
          LikelihoodEstimator* new_likelihood_estimator = new ExactLikelihoodEstimator(rng_in,
                                                                                       seed_in,
                                                                                       data_in,
                                                                                       new_factor,
                                                                                       factor_is_fixed_in_smc);
          
          full_index_vector.push_back(likelihood_estimators.size());
          without_cancelled_index_vector.push_back(likelihood_estimators.size());
          
          if (any_annealing)
          {
            likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(rng_in,
                                                                            seed_in,
                                                                            data_in,
                                                                            new_likelihood_estimator,
                                                                            power,
                                                                            annealing_variable,
                                                                            false));
          }
          else
          {
            likelihood_estimators.push_back(new_likelihood_estimator);
          }
          
        }
        else if ( current_factor.containsElementNamed("sbi_likelihood") )
        {
          SimulateModelPtr simulate_data_model_in;
          if (current_factor.containsElementNamed("simulate_data_model"))
          {
            SEXP simulate_data_model_SEXP = current_factor["simulate_data_model"];
            simulate_data_model_in = load_simulate_data_model(simulate_data_model_SEXP);
          }
          else
          {
            stop("sbi_likelihood factor must also contain simulate_data_model.");
          }
          
          LikelihoodEstimator* new_likelihood_estimator;
          
          List current_sbi = current_factor["sbi_likelihood"];
          if ( current_sbi.containsElementNamed("type") && current_sbi.containsElementNamed("variables") )
          {
            std::string type = current_sbi["type"];
            
            if ( (type=="euclidean_uniform_abc") || (type=="lp_uniform_abc") )
            {
              std::vector<std::string> data_variables_in;
              size_t number_of_points;
              std::string tolerance_variable;
              double tolerance;
              double p;
              bool parallel;
              size_t grain_size;
              
              if (type=="euclidean_uniform_abc")
              {
                List info = get_abc_euclidean_uniform_parameter_info(model_parameters,
                                                                     current_sbi,
                                                                     type);
                
                std::vector<std::string> temp_data_variables = info[0];
                data_variables_in = temp_data_variables;
                number_of_points = info[1];
                std::string temp_tolerance_variable = info[2];
                tolerance_variable = temp_tolerance_variable;
                tolerance = info[3];
                p = 2.0;
                parallel = info[4];
                grain_size = info[5];
              }
              else
              {
                List info = get_abc_lp_uniform_parameter_info(model_parameters,
                                                              current_sbi,
                                                              type);
                
                std::vector<std::string> temp_data_variables = info[0];
                data_variables_in = temp_data_variables;
                number_of_points = info[1];
                std::string temp_tolerance_variable = info[2];
                tolerance_variable = temp_tolerance_variable;
                tolerance = info[3];
                p = info[4];
                parallel = info[5];
                grain_size = info[6];
              }
              
              bool adaptive = false;
              for (auto k=sequencer_variables.begin();
                   k!=sequencer_variables.end();
                   ++k)
              {
                if (*k==tolerance_variable)
                {
                  adaptive = true;
                }
              }
              
              //std::vector<std::string> data_variables_in;
              //data_variables_in.push_back(data_variable);
              
              if (current_factor.containsElementNamed("summary_statistics"))
              {
                SEXP summary_statistics_SEXP = current_factor["summary_statistics"];
                SummaryStatisticsPtr summary_stats = load_summary_statistics(summary_statistics_SEXP);
                
                data_created_in_get_likelihood_estimators.push_back(summary_stats(*data_in));
                
                if (adaptive==false)
                {
                  new_likelihood_estimator = make_fixed_epsilon_lp_uniform_abc_likelihood(rng_in,
                                                                                          seed_in,
                                                                                          &data_created_in_get_likelihood_estimators.back(),
                                                                                          p,
                                                                                          data_variables_in,
                                                                                          "",
                                                                                          tolerance_variable,
                                                                                          tolerance,
                                                                                          simulate_data_model_in,
                                                                                          std::make_shared<Transform>(summary_stats),
                                                                                          number_of_points,
                                                                                          parallel,
                                                                                          grain_size);
                }
                else
                {
                  new_likelihood_estimator = make_varying_epsilon_lp_uniform_abc_likelihood(rng_in,
                                                                                            seed_in,
                                                                                            &data_created_in_get_likelihood_estimators.back(),
                                                                                            p,
                                                                                            data_variables_in,
                                                                                            "",
                                                                                            tolerance_variable,
                                                                                            simulate_data_model_in,
                                                                                            std::make_shared<Transform>(summary_stats),
                                                                                            number_of_points,
                                                                                            parallel,
                                                                                            grain_size);
                }
              }
              else
              {
                if (adaptive==false)
                {
                  new_likelihood_estimator = make_fixed_epsilon_lp_uniform_abc_likelihood(rng_in,
                                                                                          seed_in,
                                                                                          data_in,
                                                                                          p,
                                                                                          data_variables_in,
                                                                                          "",
                                                                                          tolerance_variable,
                                                                                          tolerance,
                                                                                          simulate_data_model_in,
                                                                                          number_of_points,
                                                                                          parallel,
                                                                                          grain_size);
                }
                else
                {
                  new_likelihood_estimator = make_varying_epsilon_lp_uniform_abc_likelihood(rng_in,
                                                                                            seed_in,
                                                                                            data_in,
                                                                                            p,
                                                                                            data_variables_in,
                                                                                            "",
                                                                                            tolerance_variable,
                                                                                            simulate_data_model_in,
                                                                                            number_of_points,
                                                                                            parallel,
                                                                                            grain_size);
                }
              }
              
            }
            else if (type=="gaussian_abc")
            {
              std::vector<std::string> data_variables_in;
              size_t number_of_points;
              std::string tolerance_variable;
              double tolerance;
              bool parallel;
              size_t grain_size;
              
              // same function as getting a Gaussian
              List info = get_abc_euclidean_uniform_parameter_info(model_parameters,
                                                                   current_sbi,
                                                                   type);
              
              std::vector<std::string> temp_data_variables = info[0];
              data_variables_in = temp_data_variables;
              number_of_points = info[1];
              std::string temp_tolerance_variable = info[2];
              tolerance_variable = temp_tolerance_variable;
              tolerance = info[3];
              parallel = info[4];
              grain_size = info[5];
              
              bool adaptive = false;
              for (auto k=sequencer_variables.begin();
                   k!=sequencer_variables.end();
                   ++k)
              {
                if (*k==tolerance_variable)
                {
                  adaptive = true;
                }
              }
              
              //std::vector<std::string> data_variables_in;
              //data_variables_in.push_back(data_variable);
              
              if (current_factor.containsElementNamed("summary_statistics"))
              {
                SEXP summary_statistics_SEXP = current_factor["summary_statistics"];
                SummaryStatisticsPtr summary_stats = load_summary_statistics(summary_statistics_SEXP);
                
                data_created_in_get_likelihood_estimators.push_back(summary_stats(*data_in));
                
                if (adaptive==false)
                {
                  new_likelihood_estimator = make_fixed_epsilon_gaussian_abc_likelihood(rng_in,
                                                                                        seed_in,
                                                                                        &data_created_in_get_likelihood_estimators.back(),
                                                                                        data_variables_in,
                                                                                        "",
                                                                                        tolerance_variable,
                                                                                        tolerance,
                                                                                        simulate_data_model_in,
                                                                                        std::make_shared<Transform>(summary_stats),
                                                                                        number_of_points,
                                                                                        parallel,
                                                                                        grain_size);
                }
                else
                {
                  new_likelihood_estimator = make_varying_epsilon_gaussian_abc_likelihood(rng_in,
                                                                                          seed_in,
                                                                                          &data_created_in_get_likelihood_estimators.back(),
                                                                                          data_variables_in,
                                                                                          "",
                                                                                          tolerance_variable,
                                                                                          simulate_data_model_in,
                                                                                          std::make_shared<Transform>(summary_stats),
                                                                                          number_of_points,
                                                                                          parallel,
                                                                                          grain_size);
                }
              }
              else
              {
                if (adaptive==false)
                {
                  new_likelihood_estimator = make_fixed_epsilon_gaussian_abc_likelihood(rng_in,
                                                                                        seed_in,
                                                                                        data_in,
                                                                                        data_variables_in,
                                                                                        "",
                                                                                        tolerance_variable,
                                                                                        tolerance,
                                                                                        simulate_data_model_in,
                                                                                        number_of_points,
                                                                                        parallel,
                                                                                        grain_size);
                }
                else
                {
                  new_likelihood_estimator = make_varying_epsilon_gaussian_abc_likelihood(rng_in,
                                                                                          seed_in,
                                                                                          data_in,
                                                                                          data_variables_in,
                                                                                          "",
                                                                                          tolerance_variable,
                                                                                          simulate_data_model_in,
                                                                                          number_of_points,
                                                                                          parallel,
                                                                                          grain_size);
                }
              }
              
            }
            else if (type=="enki_abc")
            {
              std::vector<std::string> data_variables_in;
              size_t number_of_points;
              std::string tolerance_variable;
              std::vector<double> schedule;
              size_t enki_lag;
              std::string shifter_name;
              double enki_annealing_desired_cess;
              size_t enki_number_of_bisections;
              bool enki_on_summary;
              double significance_level;
              bool parallel;
              size_t grain_size;
              
              List info = get_abc_enki_parameter_info(model_parameters,
                                                      current_sbi,
                                                      type);
              
              std::vector<std::string> temp_data_variables = info[0];
              data_variables_in = temp_data_variables;
              number_of_points = info[1];
              std::string temp_tolerance_variable = info[2];
              tolerance_variable = temp_tolerance_variable;
              std::vector<double> temp_schedule = info[3];
              schedule = temp_schedule;
              enki_lag = info[4];
              std::string temp_shifter_name = info[5];
              shifter_name = temp_shifter_name;
              enki_annealing_desired_cess = info[6];
              enki_number_of_bisections = info[7];
              enki_on_summary = info[8];
              significance_level = info[9];
              parallel = info[10];
              grain_size = info[11];
              
              EnsembleShifter* shifter = new StochasticEnsembleShifter();
              
              bool adaptive = false;
              for (auto k=sequencer_variables.begin();
                   k!=sequencer_variables.end();
                   ++k)
              {
                if (*k==tolerance_variable)
                {
                  adaptive = true;
                }
              }
              
              //std::vector<std::string> data_variables_in;
              //data_variables_in.push_back(data_variable);
              
              if (current_factor.containsElementNamed("summary_statistics"))
              {
                SEXP summary_statistics_SEXP = current_factor["summary_statistics"];
                SummaryStatisticsPtr summary_stats = load_summary_statistics(summary_statistics_SEXP);
                
                data_created_in_get_likelihood_estimators.push_back(summary_stats(*data_in));
                
                if (adaptive==false)
                {
                  new_likelihood_estimator = make_fixed_epsilon_enki_abc_likelihood(rng_in,
                                                                                    seed_in,
                                                                                    &data_created_in_get_likelihood_estimators.back(),
                                                                                    enki_lag,
                                                                                    shifter,
                                                                                    enki_annealing_desired_cess,
                                                                                    enki_number_of_bisections,
                                                                                    "",
                                                                                    tolerance_variable,
                                                                                    schedule,
                                                                                    data_variables_in,
                                                                                    simulate_data_model_in,
                                                                                    std::make_shared<Transform>(summary_stats),
                                                                                    enki_on_summary,
                                                                                    number_of_points,
                                                                                    significance_level,
                                                                                    parallel,
                                                                                    grain_size);
                }
                else
                {
                  stop("Adaptive approach to choosing epsilon not yet implemented in EnKI-ABC.");
                }
              }
              else
              {
                if (adaptive==false)
                {
                  new_likelihood_estimator = make_fixed_epsilon_enki_abc_likelihood(rng_in,
                                                                                    seed_in,
                                                                                    data_in,
                                                                                    enki_lag,
                                                                                    shifter,
                                                                                    enki_annealing_desired_cess,
                                                                                    enki_number_of_bisections,
                                                                                    "",
                                                                                    tolerance_variable,
                                                                                    schedule,
                                                                                    data_variables_in,
                                                                                    simulate_data_model_in,
                                                                                    number_of_points,
                                                                                    significance_level,
                                                                                    parallel,
                                                                                    grain_size);
                }
                else
                {
                  stop("Adaptive approach to choosing epsilon not yet implemented in EnKI-ABC.");
                }
              }
              
            }
            else
            {
              Rcout << "SBI type " << type;
              stop("SBI type unknown.");
            }
            
          }
          else
          {
            stop("Missing information for SBI in model file.");
          }
          
          full_index_vector.push_back(likelihood_estimators.size());
          without_cancelled_index_vector.push_back(likelihood_estimators.size());
          
          if (any_annealing)
          {
            likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(rng_in,
                                                                            seed_in,
                                                                            data_in,
                                                                            new_likelihood_estimator,
                                                                            power,
                                                                            annealing_variable,
                                                                            false));
          }
          else
          {
            likelihood_estimators.push_back(new_likelihood_estimator);
          }
          
        }
        else if ( current_factor.containsElementNamed("linear_gaussian_data_covariance") )
        {
          if (!current_factor.containsElementNamed("linear_gaussian_data_matrix"))
          {
            stop("linear_gaussian_data must contain both matrix and covariance.");
          }
          
          if (!current_factor.containsElementNamed("type"))
          {
            stop("linear_gaussian_data must contain a type.");
          }
          
          if (!current_factor.containsElementNamed("linear_gaussian_data_variable"))
          {
            stop("linear_gaussian_data must contain linear_gaussian_data_variable.");
          }
          
          if (!current_factor.containsElementNamed("linear_gaussian_data_state_variable"))
          {
            stop("linear_gaussian_data must contain linear_gaussian_data_state_variable.");
          }
          
          SEXP linear_gaussian_data_matrix_SEXP = current_factor["linear_gaussian_data_matrix"];
          SEXP linear_gaussian_data_covariance_SEXP = current_factor["linear_gaussian_data_covariance"];
          SEXP linear_gaussian_data_variable_SEXP = current_factor["linear_gaussian_data_variable"];
          SEXP linear_gaussian_data_state_variable_SEXP = current_factor["linear_gaussian_data_state_variable"];
          size_t type = current_factor["type"];
          ProposalKernel* proposal;
          if (type==1)
          {
            proposal = new LinearGaussianNoiseProposalKernel(load_string(linear_gaussian_data_variable_SEXP),
                                                             load_string(linear_gaussian_data_state_variable_SEXP),
                                                             load_matrix(linear_gaussian_data_matrix_SEXP),
                                                             load_matrix(linear_gaussian_data_covariance_SEXP));
          }
          else if (type==2)
          {
            proposal = new LinearGaussianNoiseFunctionProposalKernel(load_string(linear_gaussian_data_variable_SEXP),
                                                                     load_string(linear_gaussian_data_state_variable_SEXP),
                                                                     load_matrix_function(linear_gaussian_data_matrix_SEXP),
                                                                     load_matrix_function(linear_gaussian_data_covariance_SEXP));
          }
          else
          {
            stop("linear_gaussian_data specificiation is invalid.");
          }
          
          LikelihoodEstimator* new_likelihood_estimator = new ExactLikelihoodEstimator(rng_in,
                                                                                       seed_in,
                                                                                       data_in,
                                                                                       proposal,
                                                                                       factor_is_fixed_in_smc);
          
          full_index_vector.push_back(likelihood_estimators.size());
          without_cancelled_index_vector.push_back(likelihood_estimators.size());
          
          if (any_annealing)
          {
            likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(rng_in,
                                                                            seed_in,
                                                                            data_in,
                                                                            new_likelihood_estimator,
                                                                            power,
                                                                            annealing_variable,
                                                                            false));
          }
          else
          {
            likelihood_estimators.push_back(new_likelihood_estimator);
          }
        }
        else if ( current_factor.containsElementNamed("nonlinear_gaussian_data_covariance") )
        {
          if (!current_factor.containsElementNamed("nonlinear_gaussian_data_function"))
          {
            stop("nonlinear_gaussian_data must contain both function and covariance.");
          }
          
          if (!current_factor.containsElementNamed("type"))
          {
            stop("nonlinear_gaussian_data must contain a type.");
          }
          
          if (!current_factor.containsElementNamed("nonlinear_gaussian_data_variable"))
          {
            stop("nonlinear_gaussian_data must contain nonlinear_gaussian_data_variable.");
          }
          
          SEXP nonlinear_gaussian_data_function_SEXP = current_factor["nonlinear_gaussian_data_function"];
          SEXP nonlinear_gaussian_data_covariance_SEXP = current_factor["nonlinear_gaussian_data_covariance"];
          SEXP nonlinear_gaussian_data_variable_SEXP = current_factor["nonlinear_gaussian_data_variable"];
          size_t type = current_factor["type"];
          ProposalKernel* proposal;
          if (type==1)
          {
            proposal = new NonLinearGaussianNoiseProposalKernel(load_string(nonlinear_gaussian_data_variable_SEXP),
                                                                std::make_shared<Transform>(load_transform(nonlinear_gaussian_data_function_SEXP)),
                                                                load_matrix(nonlinear_gaussian_data_covariance_SEXP));
          }
          else if (type==2)
          {
            proposal = new NonLinearGaussianNoiseFunctionProposalKernel(load_string(nonlinear_gaussian_data_variable_SEXP),
                                                                        std::make_shared<Transform>(load_transform(nonlinear_gaussian_data_function_SEXP)),
                                                                        load_matrix_function(nonlinear_gaussian_data_covariance_SEXP));
          }
          else
          {
            stop("nonlinear_gaussian_data specificiation is invalid.");
          }
          
          LikelihoodEstimator* new_likelihood_estimator = new ExactLikelihoodEstimator(rng_in,
                                                                                       seed_in,
                                                                                       data_in,
                                                                                       proposal,
                                                                                       factor_is_fixed_in_smc);
          
          full_index_vector.push_back(likelihood_estimators.size());
          without_cancelled_index_vector.push_back(likelihood_estimators.size());
          
          if (any_annealing)
          {
            likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(rng_in,
                                                                            seed_in,
                                                                            data_in,
                                                                            new_likelihood_estimator,
                                                                            power,
                                                                            annealing_variable,
                                                                            false));
          }
          else
          {
            likelihood_estimators.push_back(new_likelihood_estimator);
          }
        }
        else
        {
          stop("Invalid factor.");
        }
        
      }
      else
      {
        stop("Error in factors section of model file.");
      }
      
    }
    
  }

  /*
  likelihood_estimators.push_back(new ExactLikelihoodEstimator(rng_in,
                                                               seed_in,
                                                               data_in,
                                                               prior_factors,
                                                               exact_likelihood_factors,
                                                               true));
  */
  
  
  if ( (include_priors) && (any_annealing) )
  {
    full_index_vector.push_back(likelihood_estimators.size());
    without_cancelled_index_vector.push_back(likelihood_estimators.size());
    
    ExactLikelihoodEstimator* proposal_for_evaluation = new ExactLikelihoodEstimator(rng_in,
                                                                                     seed_in,
                                                                                     data_in,
                                                                                     proposal_in,
                                                                                     true);
    likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(rng_in,
                                                                    seed_in,
                                                                    data_in,
                                                                    proposal_for_evaluation,
                                                                    second_power,
                                                                    annealing_variable,
                                                                    false));
    
  }
  
  
  full_index = new VectorSingleIndex(full_index_vector);
  without_cancelled_index = new VectorSingleIndex(without_cancelled_index_vector);
  
  /*
  Rcout << full_index_vector.size() << std::endl;
  Rcout << without_cancelled_index_vector.size() << std::endl;
  Rcout << likelihood_estimators.size() << std::endl;
  */
  
  return likelihood_estimators;
  
  //return std::vector<LikelihoodEstimator*>();
}


KalmanUpdater* get_kalman_updater(const List &model,
                                  const List &model_parameters)
{
  
  if ( model.containsElementNamed("factor") )
  {
    List factors = model["factor"];
    
    for (size_t i=0; i<factors.size(); ++i)
    {
      
      if (Rf_isNewList(factors[i]))
      {
        List current_factor = factors[i];
        
        if ( current_factor.containsElementNamed("linear_gaussian_data_covariance") )
        {
          if (!current_factor.containsElementNamed("linear_gaussian_data_matrix"))
          {
            stop("linear_gaussian must contain both matrix and covariance.");
          }
          
          if (!current_factor.containsElementNamed("type"))
          {
            stop("linear_gaussian must contain a type.");
          }
          
          //if (!current_factor.containsElementNamed("linear_gaussian_data_variable"))
          //{
          //  stop("linear_gaussian_data must contain linear_gaussian_data_variable.");
          //}
          
          //if (!current_factor.containsElementNamed("linear_gaussian_data_state_variable"))
          //{
          //  stop("linear_gaussian_data must contain linear_gaussian_data_state_variable.");
          //}
          
          SEXP linear_gaussian_data_matrix_SEXP = current_factor["linear_gaussian_data_matrix"];
          SEXP linear_gaussian_data_covariance_SEXP = current_factor["linear_gaussian_data_covariance"];
          //SEXP linear_gaussian_data_variable_SEXP = current_factor["linear_gaussian_data_variable"];
          //SEXP linear_gaussian_data_state_variable_SEXP = current_factor["linear_gaussian_data_state_variable"];
          size_t type = current_factor["type"];

          if (type==1)
          {
            return new ExactKalmanUpdater(load_matrix(linear_gaussian_data_matrix_SEXP),
                                          load_matrix(linear_gaussian_data_covariance_SEXP));
          }
          else if (type==2)
          {
            return new ExactKalmanUpdater(load_matrix_function(linear_gaussian_data_matrix_SEXP),
                                          load_matrix_function(linear_gaussian_data_covariance_SEXP));
          }
          else
          {
            stop("linear_gaussian specificiation is invalid.");
          }
          
        }
        else if ( current_factor.containsElementNamed("nonlinear_gaussian_data_covariance") )
        {
          stop("Unscented Kalman filter not yet supported.");
          
          /*
           if (!current_factor.containsElementNamed("nonlinear_gaussian_data_function"))
           {
           stop("nonlinear_gaussian_data must contain both function and covariance.");
           }
           
           if (!current_factor.containsElementNamed("type"))
           {
           stop("nonlinear_gaussian_data must contain a type.");
           }
           
           if (!current_factor.containsElementNamed("nonlinear_gaussian_data_variable"))
           {
           stop("nonlinear_gaussian_data must contain nonlinear_gaussian_data_variable.");
           }
           
           SEXP nonlinear_gaussian_data_function_SEXP = current_factor["nonlinear_gaussian_data_function"];
           SEXP nonlinear_gaussian_data_covariance_SEXP = current_factor["nonlinear_gaussian_data_covariance"];
           SEXP nonlinear_gaussian_data_variable_SEXP = current_factor["nonlinear_gaussian_data_variable"];
           size_t type = current_factor["type"];
           ProposalKernel* proposal;
           if (type==1)
           {
           proposal = new NonLinearGaussianNoiseProposalKernel(load_string(nonlinear_gaussian_data_variable_SEXP),
           std::make_shared<Transform>(load_transform(nonlinear_gaussian_data_function_SEXP)),
           load_matrix(nonlinear_gaussian_data_covariance_SEXP));
           }
           else if (type==2)
           {
           proposal = new NonLinearGaussianNoiseFunctionProposalKernel(load_string(nonlinear_gaussian_data_variable_SEXP),
           std::make_shared<Transform>(load_transform(nonlinear_gaussian_data_function_SEXP)),
           load_matrix_function(nonlinear_gaussian_data_covariance_SEXP));
           }
           else
           {
           stop("nonlinear_gaussian_data specificiation is invalid.");
           }
           */
        }
        
      }
      else
      {
        stop("Error in factors section of model file.");
      }
      
    }
    
  }
  else
  {
    stop("No factors found in model file.");
  }

  stop("Valid Kalman updater not found.");
}

KalmanPredictor* get_kalman_predictor(const List &model,
                                      const List &model_parameters)
{
  
  if ( model.containsElementNamed("transition_model") )
  {
    List factors = model["transition_model"];
    
    for (size_t i=0; i<factors.size(); ++i)
    {
      
      if (Rf_isNewList(factors[i]))
      {
        List current_factor = factors[i];
        
        if ( current_factor.containsElementNamed("linear_gaussian_transition_covariance") )
        {
          if (!current_factor.containsElementNamed("linear_gaussian_transition_matrix"))
          {
            stop("linear_gaussian must contain both matrix and covariance.");
          }
          
          if (!current_factor.containsElementNamed("type"))
          {
            stop("linear_gaussian must contain a type.");
          }
          
          //if (!current_factor.containsElementNamed("linear_gaussian_data_variable"))
          //{
          //  stop("linear_gaussian_data must contain linear_gaussian_data_variable.");
          //}
          
          //if (!current_factor.containsElementNamed("linear_gaussian_data_state_variable"))
          //{
          //  stop("linear_gaussian_data must contain linear_gaussian_data_state_variable.");
          //}
          
          SEXP linear_gaussian_transition_matrix_SEXP = current_factor["linear_gaussian_transition_matrix"];
          SEXP linear_gaussian_transition_covariance_SEXP = current_factor["linear_gaussian_transition_covariance"];
          //SEXP linear_gaussian_data_variable_SEXP = current_factor["linear_gaussian_data_variable"];
          //SEXP linear_gaussian_data_state_variable_SEXP = current_factor["linear_gaussian_data_state_variable"];
          size_t type = current_factor["type"];
          
          if (type==1)
          {
            return new ExactKalmanPredictor(load_matrix(linear_gaussian_transition_matrix_SEXP),
                                            load_matrix(linear_gaussian_transition_covariance_SEXP));
          }
          else if (type==2)
          {
            return new ExactKalmanPredictor(load_matrix_function(linear_gaussian_transition_matrix_SEXP),
                                            load_matrix_function(linear_gaussian_transition_covariance_SEXP));
          }
          else
          {
            stop("linear_gaussian specificiation is invalid.");
          }
          
        }
        else if ( current_factor.containsElementNamed("nonlinear_gaussian_transition_covariance") )
        {
          stop("Unscented Kalman filter not yet supported.");
          
          /*
           if (!current_factor.containsElementNamed("nonlinear_gaussian_data_function"))
           {
           stop("nonlinear_gaussian_data must contain both function and covariance.");
           }
           
           if (!current_factor.containsElementNamed("type"))
           {
           stop("nonlinear_gaussian_data must contain a type.");
           }
           
           if (!current_factor.containsElementNamed("nonlinear_gaussian_data_variable"))
           {
           stop("nonlinear_gaussian_data must contain nonlinear_gaussian_data_variable.");
           }
           
           SEXP nonlinear_gaussian_data_function_SEXP = current_factor["nonlinear_gaussian_data_function"];
           SEXP nonlinear_gaussian_data_covariance_SEXP = current_factor["nonlinear_gaussian_data_covariance"];
           SEXP nonlinear_gaussian_data_variable_SEXP = current_factor["nonlinear_gaussian_data_variable"];
           size_t type = current_factor["type"];
           ProposalKernel* proposal;
           if (type==1)
           {
           proposal = new NonLinearGaussianNoiseProposalKernel(load_string(nonlinear_gaussian_data_variable_SEXP),
           std::make_shared<Transform>(load_transform(nonlinear_gaussian_data_function_SEXP)),
           load_matrix(nonlinear_gaussian_data_covariance_SEXP));
           }
           else if (type==2)
           {
           proposal = new NonLinearGaussianNoiseFunctionProposalKernel(load_string(nonlinear_gaussian_data_variable_SEXP),
           std::make_shared<Transform>(load_transform(nonlinear_gaussian_data_function_SEXP)),
           load_matrix_function(nonlinear_gaussian_data_covariance_SEXP));
           }
           else
           {
           stop("nonlinear_gaussian_data specificiation is invalid.");
           }
           */
        }
        
      }
      else
      {
        stop("Error in transition_models section of model file.");
      }
      
    }
    
  }
  else
  {
    stop("No transition_models found in model file.");
  }
  
  stop("Valid Kalman predictor not found.");
}

List get_prior_mean_and_covariance(const List &model,
                                   const List &model_parameters)
{
  if ( model.containsElementNamed("factor") )
  {
    
    List factors = model["factor"];
    
    for (size_t i=0; i<factors.size(); ++i)
    {
      
      if (Rf_isNewList(factors[i]))
      {
        
        List current_factor = factors[i];
        if ( current_factor.containsElementNamed("prior") )
        {
          DistributionFactor* new_factor;
          
          if (Rf_isNewList(current_factor["prior"]))
          {
            
            List current_prior = current_factor["prior"];
            if ( current_prior.containsElementNamed("type") && current_prior.containsElementNamed("variables") )
            {
              std::string type = current_prior["type"];
              
              if (type=="norm")
              {
                List info = get_single_variable_two_double_parameter_info(model_parameters,
                                                                          current_prior,
                                                                          type);
                
                std::string variable = info[0];
                double mean = info[1];
                double sd = info[2];
                
                arma::colvec mean_vec(1);
                mean_vec[0] = mean;
                
                arma::colvec cov(1,1);
                cov(0,0) = sd*sd;
                
                return List::create(mean_vec,cov);

                
              }
              else if (type=="mvnorm")
              {
                List info = get_single_variable_vector_and_matrix_parameter_info(model_parameters,
                                                                                 current_prior,
                                                                                 type);
                
                std::string variable = info[0];
                arma::colvec mean = info[1];
                arma::mat cov = info[2];
                
                return List::create(mean,cov);
                
              }
              
            }
            else
            {
              stop("Missing information for prior in model file.");
            }
            
          }
          else
          {
            stop("Error in factors section of model file.");
          }
          
        }
        else
        {
          stop("Invalid factor.");
        }
        
      }
      else
      {
        stop("Error in factors section of model file.");
      }
      
    }
    
  }
  else
  {
    stop("No factors found in model.");
  }
  
  stop("Valid prior mean and covariance not found.");
}



std::vector<MeasurementCovarianceEstimator*> get_measurement_covariance_estimators(RandomNumberGenerator* rng_in,
                                                                                   size_t* seed_in,
                                                                                   Data* data_in,
                                                                                   const List &model,
                                                                                   const List &model_parameters,
                                                                                   const std::vector<std::string> &sequencer_types,
                                                                                   const std::vector<std::string> &sequencer_variables,
                                                                                   const std::vector<std::vector<double>>  &sequencer_schedules,
                                                                                   const List &enk_likelihood_index_method,
                                                                                   std::shared_ptr<Transform> transform,
                                                                                   std::vector<Data> &data_created_in_get_measurement_covariance_estimators)
{
  std::vector<size_t> likelihood_indices;
  if (enk_likelihood_index_method.containsElementNamed("func"))
  {
    NumericVector enk_likelihood_indices;
    SEXP enk_likelihood_index_SEXP = enk_likelihood_index_method["func"];
    enk_likelihood_indices = load_numeric_vector(enk_likelihood_index_SEXP);
    
    likelihood_indices.reserve(enk_likelihood_indices.length());
    for (size_t i=0; i<enk_likelihood_indices.length(); ++i)
    {
      likelihood_indices.push_back(size_t(enk_likelihood_indices[i]));
    }
  }

  data_created_in_get_measurement_covariance_estimators.clear();
  
  std::vector<MeasurementCovarianceEstimator*> measurement_covariance_estimators;
  //std::vector<DistributionFactor*> numerator_distribution_factors;
  //std::vector<LikelihoodFactor*> numerator_likelihood_factors;
  
  //std::vector<DistributionFactor*> prior_factors;
  //std::vector<LikelihoodFactor*> exact_likelihood_factors;
  
  //std::vector<IndependentProposalKernel*> numerator_proposal_factors;
  
  if ( model.containsElementNamed("factor") )
  {
    List factors = model["factor"];
    
    std::vector<size_t> indices_to_use;
    if (likelihood_indices.size()==0)
    {
      indices_to_use.reserve(factors.length());
      for (size_t i=0; i<factors.length(); ++i)
      {
        indices_to_use.push_back(i);
      }
    }
    else
    {
      indices_to_use = likelihood_indices;
    }
    
    for (size_t i=0; i<indices_to_use.size(); ++i)
    {
      if (indices_to_use[i]>=factors.length())
      {
        stop("Error in enk_likelihood_index: index does not index a specified factor.");
      }
      
      if (Rf_isNewList(factors[indices_to_use[i]]))
      {
        List current_factor = factors[indices_to_use[i]];
        SimulateModelPtr simulate_data_model_in;
        if (current_factor.containsElementNamed("simulate_data_model"))
        {
          SEXP simulate_data_model_SEXP = current_factor["simulate_data_model"];
          simulate_data_model_in = load_simulate_data_model(simulate_data_model_SEXP);
          
          MeasurementCovarianceEstimator* new_measurement_covariance_estimator;
          
          if (current_factor.containsElementNamed("summary_statistics"))
          {
            SEXP summary_statistics_SEXP = current_factor["summary_statistics"];
            SummaryStatisticsPtr summary_stats = load_summary_statistics(summary_statistics_SEXP);
            
            data_created_in_get_measurement_covariance_estimators.push_back(summary_stats(*data_in));
            
            new_measurement_covariance_estimator = new GenericMeasurementCovarianceEstimator(rng_in,
                                                                                             seed_in,
                                                                                             &data_created_in_get_measurement_covariance_estimators.back(),
                                                                                             transform,
                                                                                             std::make_shared<Transform>(summary_stats),
                                                                                             simulate_data_model_in);
          }
          else
          {
            
            new_measurement_covariance_estimator = new GenericMeasurementCovarianceEstimator(rng_in,
                                                                                             seed_in,
                                                                                             data_in,
                                                                                             transform,
                                                                                             NULL,
                                                                                             simulate_data_model_in);
          }
          
          measurement_covariance_estimators.push_back(new_measurement_covariance_estimator);
        }
        else if ( current_factor.containsElementNamed("linear_gaussian_data_covariance") )
        {
          if (!current_factor.containsElementNamed("linear_gaussian_data_matrix"))
          {
            stop("linear_gaussian_data must contain both matrix and covariance.");
          }
          
          if (!current_factor.containsElementNamed("type"))
          {
            stop("linear_gaussian_data must contain a type.");
          }
          
          if (!current_factor.containsElementNamed("linear_gaussian_data_variable"))
          {
            stop("linear_gaussian_data must contain linear_gaussian_data_variable.");
          }
          
          if (!current_factor.containsElementNamed("linear_gaussian_data_state_variable"))
          {
            stop("linear_gaussian_data must contain linear_gaussian_data_state_variable.");
          }
          
          SEXP linear_gaussian_data_matrix_SEXP = current_factor["linear_gaussian_data_matrix"];
          SEXP linear_gaussian_data_covariance_SEXP = current_factor["linear_gaussian_data_covariance"];
          SEXP linear_gaussian_data_variable_SEXP = current_factor["linear_gaussian_data_variable"];
          SEXP linear_gaussian_data_state_variable_SEXP = current_factor["linear_gaussian_data_state_variable"];
          size_t type = current_factor["type"];
          
          MeasurementCovarianceEstimator* new_measurement_covariance_estimator;
          
          if (current_factor.containsElementNamed("summary_statistics"))
          {
            SEXP summary_statistics_SEXP = current_factor["summary_statistics"];
            SummaryStatisticsPtr summary_stats = load_summary_statistics(summary_statistics_SEXP);
            
            data_created_in_get_measurement_covariance_estimators.push_back(summary_stats(*data_in));
            
            if (type==1)
            {
              new_measurement_covariance_estimator = new DirectLinearGaussianMeasurementCovarianceEstimator(rng_in,
                                                                                                            seed_in,
                                                                                                            &data_created_in_get_measurement_covariance_estimators.back(),
                                                                                                            transform,
                                                                                                            std::make_shared<Transform>(summary_stats),
                                                                                                            load_matrix(linear_gaussian_data_matrix_SEXP),
                                                                                                            load_matrix(linear_gaussian_data_covariance_SEXP),
                                                                                                            load_string(linear_gaussian_data_variable_SEXP),
                                                                                                            load_string(linear_gaussian_data_state_variable_SEXP));
            }
            else if (type==2)
            {
              new_measurement_covariance_estimator = new DirectLinearGaussianMeasurementCovarianceEstimator(rng_in,
                                                                                                            seed_in,
                                                                                                            &data_created_in_get_measurement_covariance_estimators.back(),
                                                                                                            transform,
                                                                                                            std::make_shared<Transform>(summary_stats),
                                                                                                            load_matrix_function(linear_gaussian_data_matrix_SEXP),
                                                                                                            load_matrix_function(linear_gaussian_data_covariance_SEXP),
                                                                                                            load_string(linear_gaussian_data_variable_SEXP),
                                                                                                            load_string(linear_gaussian_data_state_variable_SEXP));
            }
            else
            {
              stop("linear_gaussian_data specificiation is invalid.");
            }
              
          }
          else
          {
            
            if (type==1)
            {
              new_measurement_covariance_estimator = new DirectLinearGaussianMeasurementCovarianceEstimator(rng_in,
                                                                                                            seed_in,
                                                                                                            data_in,
                                                                                                            transform,
                                                                                                            NULL,
                                                                                                            load_matrix(linear_gaussian_data_matrix_SEXP),
                                                                                                            load_matrix(linear_gaussian_data_covariance_SEXP),
                                                                                                            load_string(linear_gaussian_data_variable_SEXP),
                                                                                                            load_string(linear_gaussian_data_state_variable_SEXP));
            }
            else if (type==2)
            {
              new_measurement_covariance_estimator = new DirectLinearGaussianMeasurementCovarianceEstimator(rng_in,
                                                                                                            seed_in,
                                                                                                            data_in,
                                                                                                            transform,
                                                                                                            NULL,
                                                                                                            load_matrix_function(linear_gaussian_data_matrix_SEXP),
                                                                                                            load_matrix_function(linear_gaussian_data_covariance_SEXP),
                                                                                                            load_string(linear_gaussian_data_variable_SEXP),
                                                                                                            load_string(linear_gaussian_data_state_variable_SEXP));
            }
            else
            {
              stop("linear_gaussian_data specificiation is invalid.");
            }
            
            
          }
          
          measurement_covariance_estimators.push_back(new_measurement_covariance_estimator);
        }
        else if ( current_factor.containsElementNamed("nonlinear_gaussian_data_covariance") )
        {
          if (!current_factor.containsElementNamed("nonlinear_gaussian_data_function"))
          {
            stop("nonlinear_gaussian_data must contain both function and covariance.");
          }
          
          if (!current_factor.containsElementNamed("type"))
          {
            stop("nonlinear_gaussian_data must contain a type.");
          }
          
          if (!current_factor.containsElementNamed("nonlinear_gaussian_data_variable"))
          {
            stop("nonlinear_gaussian_data must contain nonlinear_gaussian_data_variable.");
          }
          
          SEXP nonlinear_gaussian_data_function_SEXP = current_factor["nonlinear_gaussian_data_function"];
          SEXP nonlinear_gaussian_data_covariance_SEXP = current_factor["nonlinear_gaussian_data_covariance"];
          SEXP nonlinear_gaussian_data_variable_SEXP = current_factor["nonlinear_gaussian_data_variable"];
          size_t type = current_factor["type"];
          
          TransformPtr nonlinear_gaussian_data_function = load_transform(nonlinear_gaussian_data_function_SEXP);
          
          MeasurementCovarianceEstimator* new_measurement_covariance_estimator;
          
          if (current_factor.containsElementNamed("summary_statistics"))
          {
            SEXP summary_statistics_SEXP = current_factor["summary_statistics"];
            SummaryStatisticsPtr summary_stats = load_summary_statistics(summary_statistics_SEXP);
            
            data_created_in_get_measurement_covariance_estimators.push_back(summary_stats(*data_in));
            
            if (type==1)
            {
              new_measurement_covariance_estimator = new DirectNonLinearGaussianMeasurementCovarianceEstimator(rng_in,
                                                                                                              seed_in,
                                                                                                              &data_created_in_get_measurement_covariance_estimators.back(),
                                                                                                               transform,
                                                                                                              std::make_shared<Transform>(summary_stats),
                                                                                                               std::make_shared<Transform>(nonlinear_gaussian_data_function),
                                                                                                               load_matrix(nonlinear_gaussian_data_covariance_SEXP),
                                                                                                               load_string(nonlinear_gaussian_data_variable_SEXP));
            }
            else if (type==2)
            {
              new_measurement_covariance_estimator = new DirectNonLinearGaussianMeasurementCovarianceEstimator(rng_in,
                                                                                                               seed_in,
                                                                                                               &data_created_in_get_measurement_covariance_estimators.back(),
                                                                                                               transform,
                                                                                                               std::make_shared<Transform>(summary_stats),
                                                                                                               std::make_shared<Transform>(nonlinear_gaussian_data_function),
                                                                                                               load_matrix_function(nonlinear_gaussian_data_covariance_SEXP),
                                                                                                               load_string(nonlinear_gaussian_data_variable_SEXP));
            }
            else
            {
              stop("nonlinear_gaussian_data specificiation is invalid.");
            }
            
          }
          else
          {
            
            if (type==1)
            {
              new_measurement_covariance_estimator = new DirectNonLinearGaussianMeasurementCovarianceEstimator(rng_in,
                                                                                                               seed_in,
                                                                                                               data_in,
                                                                                                               transform,
                                                                                                               NULL,
                                                                                                               std::make_shared<Transform>(nonlinear_gaussian_data_function),
                                                                                                               load_matrix(nonlinear_gaussian_data_covariance_SEXP),
                                                                                                               load_string(nonlinear_gaussian_data_variable_SEXP));
            }
            else if (type==2)
            {
              new_measurement_covariance_estimator = new DirectNonLinearGaussianMeasurementCovarianceEstimator(rng_in,
                                                                                                               seed_in,
                                                                                                               data_in,
                                                                                                               transform,
                                                                                                               NULL,
                                                                                                               std::make_shared<Transform>(nonlinear_gaussian_data_function),
                                                                                                               load_matrix_function(nonlinear_gaussian_data_covariance_SEXP),
                                                                                                               load_string(nonlinear_gaussian_data_variable_SEXP));
            }
            else
            {
              stop("nonlinear_gaussian_data specificiation is invalid.");
            }
            
            
          }
          
          measurement_covariance_estimators.push_back(new_measurement_covariance_estimator);
        }
        
      }
      //else
      //{
      //  stop("For use in an ensemble Kalman method, the factor must be either linear Gaussian, non-linear Gaussian, or simulate_model.");
      //}
      
    }
    
  }
  
   /*
   likelihood_estimators.push_back(new ExactLikelihoodEstimator(rng_in,
   seed_in,
   data_in,
   prior_factors,
   exact_likelihood_factors,
   true));
   */
  
   /*
   Rcout << full_index_vector.size() << std::endl;
   Rcout << without_cancelled_index_vector.size() << std::endl;
   Rcout << likelihood_estimators.size() << std::endl;
   */
  
  return measurement_covariance_estimators;
  
  //return std::vector<LikelihoodEstimator*>();
}

IndependentProposalKernel* get_proposal(const List &model,
                                        const List &model_parameters,
                                        const Data* data)
{
  IndependentProposalKernel* proposal = NULL;
  std::vector<IndependentProposalKernel*> proposals;
  
  if ( model.containsElementNamed("importance_proposal") )
  {
    List importance_proposals = model["importance_proposal"];
    
    for (size_t i=0; i<importance_proposals.size(); ++i)
    {
      if (Rf_isNewList(importance_proposals[i]))
      {
        List current_proposal = importance_proposals[i];
        if ( current_proposal.containsElementNamed("importance_proposal") )
        {
          List proposal_info = current_proposal["importance_proposal"];
          
          if (Rf_isNewList(proposal_info))
          {
            if ( proposal_info.containsElementNamed("type") && proposal_info.containsElementNamed("variables") )
            {
              std::string type = proposal_info["type"];
              if (type=="norm")
              {
                List info = get_single_variable_two_double_parameter_info(model_parameters,
                                                                          proposal_info,
                                                                          type);
                
                std::string variable = info[0];
                double mean = info[1];
                double sd = info[2];
                
                proposal = new GaussianIndependentProposalKernel(variable,
                                                                 mean,
                                                                 sd);
                
              }
              else if (type=="mvnorm")
              {
                List info = get_single_variable_vector_and_matrix_parameter_info(model_parameters,
                                                                                 proposal_info,
                                                                                 type);
                
                std::string variable = info[0];
                arma::colvec mean = info[1];
                arma::mat cov = info[2];
                
                proposal = new GaussianIndependentProposalKernel(variable,
                                                                 mean,
                                                                 cov);
                
              }
              else if (type=="lnorm")
              {
                List info = get_single_variable_two_double_parameter_info(model_parameters,
                                                                          proposal_info,
                                                                          type);
                
                std::string variable = info[0];
                double mean = info[1];
                double sd = info[2];
                
                proposal = new LogGaussianIndependentProposalKernel(variable,
                                                                    mean,
                                                                    sd);
                
              }
              else if (type=="mvlnorm")
              {
                List info = get_single_variable_vector_and_matrix_parameter_info(model_parameters,
                                                                                 proposal_info,
                                                                                 type);
                
                std::string variable = info[0];
                arma::colvec mean = info[1];
                arma::mat cov = info[2];
                
                proposal = new LogGaussianIndependentProposalKernel(variable,
                                                                    mean,
                                                                    cov);
                
              }
              else if (type=="gamma")
              {
                List info = get_single_variable_two_double_parameter_info(model_parameters,
                                                                   proposal_info,
                                                                   type);
                
                std::string variable = info[0];
                double shape = info[1];
                double rate = info[2];
                
                proposal = new GammaIndependentProposalKernel(variable,
                                                              shape,
                                                              rate);
                
              }
              else
              {
                Rcout << "Proposal type " << type;
                stop("Proposal type unknown");
              }
            }
            else
            {
              stop("Missing information for proposal in model file.");
            }
            
          }
          else
          {
            stop("Error in proposal section of model file.");
          }
        }
        else if ( current_proposal.containsElementNamed("evaluate_log_importance_proposal") && current_proposal.containsElementNamed("simulate_importance_proposal") && current_proposal.containsElementNamed("type") )
        {
          SEXP evaluate_log_importance_proposal_SEXP = current_proposal["evaluate_log_importance_proposal"];
          SEXP simulate_importance_proposal_SEXP = current_proposal["simulate_importance_proposal"];
          size_t type = current_proposal["type"];
          if (type==1)
          {
            proposal = new CustomDistributionProposalKernel(load_simulate_distribution(simulate_importance_proposal_SEXP),
                                                            load_evaluate_log_distribution(evaluate_log_importance_proposal_SEXP));
          }
          else if (type==2)
          {
            proposal = new CustomIndependentProposalKernel(load_simulate_independent_proposal(simulate_importance_proposal_SEXP),
                                                            load_evaluate_log_independent_proposal(evaluate_log_importance_proposal_SEXP));
          }
          else if (type==3)
          {
            proposal = new CustomGuidedDistributionProposalKernel(load_simulate_guided_distribution(simulate_importance_proposal_SEXP),
                                                                  load_evaluate_log_guided_distribution(evaluate_log_importance_proposal_SEXP),data);
          }
          else if (type==4)
          {
            proposal = new CustomGuidedIndependentProposalKernel(load_simulate_guided_independent_proposal(simulate_importance_proposal_SEXP),
                                                                  load_evaluate_log_guided_independent_proposal(evaluate_log_importance_proposal_SEXP),data);
          }
          else
          {
            stop("Functions read in for proposal are invalid.");
          }
        }
        else
        {
          stop("Invalid importance proposal. Maybe you specified a method to simulate from the proposal, but no method to evaluate it?");
        }
        
      }
      else
      {
        stop("Error in importance proposals section of model file.");
      }
      
      if (proposal!=NULL)
      {
        proposals.push_back(proposal);
      }
    }
  }
  
  if (proposal==NULL)
  {
    stop("No proposal specified.");
  }
  else
  {
    if (proposals.size()==1)
    {
      return proposal;
    }
    else
    {
      return new CompositeIndependentProposalKernel(proposals);
    }
  }
}

IndependentProposalKernel* get_prior_as_proposal(const List &model,
                                                 const List &model_parameters)
{
  IndependentProposalKernel* proposal = NULL;
  std::vector<IndependentProposalKernel*> proposals;
  
  if ( model.containsElementNamed("factor") )
  {
    List factors = model["factor"];
    
    for (size_t i=0; i<factors.size(); ++i)
    {
      if (Rf_isNewList(factors[i]))
      {
        List current_factor = factors[i];
        if ( current_factor.containsElementNamed("prior") )
        {
          if (Rf_isNewList(current_factor["prior"]))
          {
            List current_prior = current_factor["prior"];
            if ( current_prior.containsElementNamed("type") && current_prior.containsElementNamed("variables") )
            {
              std::string type = current_prior["type"];
              
              if (type=="norm")
              {
                List info = get_single_variable_two_double_parameter_info(model_parameters,
                                                                          current_prior,
                                                                          type);
                
                std::string variable = info[0];
                double mean = info[1];
                double sd = info[2];
                
                proposal = new GaussianIndependentProposalKernel(variable,
                                                                 mean,
                                                                 sd);
                
                if (proposal!=NULL)
                {
                  proposals.push_back(proposal);
                }
                
              }
              else if (type=="mvnorm")
              {
                List info = get_single_variable_vector_and_matrix_parameter_info(model_parameters,
                                                                                 current_prior,
                                                                                 type);
                
                std::string variable = info[0];
                arma::colvec mean = info[1];
                arma::mat cov = info[2];
                
                proposal = new GaussianIndependentProposalKernel(variable,
                                                                 mean,
                                                                 cov);
                
                if (proposal!=NULL)
                {
                  proposals.push_back(proposal);
                }
                
              }
              else if (type=="lnorm")
              {
                List info = get_single_variable_two_double_parameter_info(model_parameters,
                                                                          current_prior,
                                                                          type);
                
                std::string variable = info[0];
                double mean = info[1];
                double sd = info[2];
                
                proposal = new LogGaussianIndependentProposalKernel(variable,
                                                                    mean,
                                                                    sd);
                
                if (proposal!=NULL)
                {
                  proposals.push_back(proposal);
                }
                
              }
              else if (type=="mvlnorm")
              {
                List info = get_single_variable_vector_and_matrix_parameter_info(model_parameters,
                                                                                 current_prior,
                                                                                 type);
                
                std::string variable = info[0];
                arma::colvec mean = info[1];
                arma::mat cov = info[2];
                
                proposal = new LogGaussianIndependentProposalKernel(variable,
                                                                    mean,
                                                                    cov);
                
                if (proposal!=NULL)
                {
                  proposals.push_back(proposal);
                }
                
              }
              else if (type=="gamma")
              {
                List info = get_single_variable_two_double_parameter_info(model_parameters,
                                                                          current_prior,
                                                                          type);
                
                std::string variable = info[0];
                double shape = info[1];
                double rate = info[2];
                
                proposal = new GammaIndependentProposalKernel(variable,
                                                              shape,
                                                              rate);
                
                if (proposal!=NULL)
                {
                  proposals.push_back(proposal);
                }
                
              }
              else
              {
                Rcout << "Prior type " << type;
                stop("Prior type unknown");
              }
            }
            else
            {
              stop("Missing information for prior in model file.");
            }
          }
          else
          {
            stop("Error in factors section of model file.");
          }
        }
        else if ((current_factor.containsElementNamed("evaluate_log_prior")) && (current_factor.containsElementNamed("simulate_prior")))
        {
          SEXP evaluate_log_prior_SEXP = current_factor["evaluate_log_prior"];
          
          SEXP simulate_prior_SEXP = current_factor["simulate_prior"];
          proposal = new CustomDistributionProposalKernel(load_simulate_distribution(simulate_prior_SEXP),
                                                          load_evaluate_log_distribution(evaluate_log_prior_SEXP));
          
          if (proposal!=NULL)
          {
            proposals.push_back(proposal);
          }
        }
        
      }
      else
      {
        stop("Error in factors section of model file.");
      }
      
    }
  }
  
  if (proposal==NULL)
  {
    stop("No suitable prior to use as proposal.");
  }
  
  if (proposal==NULL)
  {
    stop("No suitable prior to use as proposal.");
  }
  else
  {
    if (proposals.size()==1)
    {
      return proposal;
    }
    else
    {
      return new CompositeIndependentProposalKernel(proposals);
    }
  }
}

IndependentProposalKernel* get_prior_as_simulate_only_proposal(const List &model,
                                                               const List &model_parameters)
{
  IndependentProposalKernel* proposal = NULL;
  std::vector<IndependentProposalKernel*> proposals;
  
  if ( model.containsElementNamed("factor") )
  {
    List factors = model["factor"];
    
    for (size_t i=0; i<factors.size(); ++i)
    {
      if (Rf_isNewList(factors[i]))
      {
        List current_factor = factors[i];
        if ( current_factor.containsElementNamed("prior") )
        {
          if (Rf_isNewList(current_factor["prior"]))
          {
            List current_prior = current_factor["prior"];
            if ( current_prior.containsElementNamed("type") && current_prior.containsElementNamed("variables") )
            {
              std::string type = current_prior["type"];
              
              if (type=="norm")
              {
                List info = get_single_variable_two_double_parameter_info(model_parameters,
                                                                   current_prior,
                                                                   type);
                
                std::string variable = info[0];
                double mean = info[1];
                double sd = info[2];
                
                proposal = new GaussianIndependentProposalKernel(variable,
                                                                 mean,
                                                                 sd);
                
                if (proposal!=NULL)
                {
                  proposals.push_back(proposal);
                }
                
              }
              else if (type=="mvnorm")
              {
                List info = get_single_variable_vector_and_matrix_parameter_info(model_parameters,
                                                                                 current_prior,
                                                                                 type);
                
                std::string variable = info[0];
                arma::colvec mean = info[1];
                arma::mat cov = info[2];
                
                proposal = new GaussianIndependentProposalKernel(variable,
                                                                 mean,
                                                                 cov);
                
                if (proposal!=NULL)
                {
                  proposals.push_back(proposal);
                }
                
              }
              else if (type=="lnorm")
              {
                List info = get_single_variable_two_double_parameter_info(model_parameters,
                                                                   current_prior,
                                                                   type);
                
                std::string variable = info[0];
                double mean = info[1];
                double sd = info[2];
                
                proposal = new LogGaussianIndependentProposalKernel(variable,
                                                                    mean,
                                                                    sd);
                
                if (proposal!=NULL)
                {
                  proposals.push_back(proposal);
                }
                
              }
              else if (type=="mvlnorm")
              {
                List info = get_single_variable_vector_and_matrix_parameter_info(model_parameters,
                                                                                 current_prior,
                                                                                 type);
                
                std::string variable = info[0];
                arma::colvec mean = info[1];
                arma::mat cov = info[2];
                
                proposal = new LogGaussianIndependentProposalKernel(variable,
                                                                    mean,
                                                                    cov);
                
                if (proposal!=NULL)
                {
                  proposals.push_back(proposal);
                }
                
              }
              else if (type=="gamma")
              {
                List info = get_single_variable_two_double_parameter_info(model_parameters,
                                                                   current_prior,
                                                                   type);
                
                std::string variable = info[0];
                double shape = info[1];
                double rate = info[2];
                
                proposal = new GammaIndependentProposalKernel(variable,
                                                              shape,
                                                              rate);
                
                if (proposal!=NULL)
                {
                  proposals.push_back(proposal);
                }
                
              }
              else
              {
                Rcout << "Prior type " << type;
                stop("Prior type unknown");
              }
            }
            else
            {
              stop("Missing information for prior in model file.");
            }
          }
          else
          {
            stop("Error in factors section of model file.");
          }
        }
        else if (current_factor.containsElementNamed("simulate_prior"))
        {
          SEXP simulate_prior_SEXP = current_factor["simulate_prior"];
          proposal = new CustomDistributionProposalKernel(load_simulate_distribution(simulate_prior_SEXP));
          
          if (proposal!=NULL)
          {
            proposals.push_back(proposal);
          }
        }
        
      }
      else
      {
        stop("Error in factors section of model file.");
      }
      
    }
  }
  
  if (proposal==NULL)
  {
    stop("No suitable prior to use as proposal.");
  }
  else
  {
    if (proposals.size()==1)
    {
      return proposal;
    }
    else
    {
      return new CompositeIndependentProposalKernel(proposals);
    }
  }
}

ProposalKernel* get_mh_proposal(const List &current_proposal,
                                const List &model_parameters,
                                const Data* data)
{
  ProposalKernel* proposal;
  if ( current_proposal.containsElementNamed("mh_proposal") )
  {
    List proposal_info = current_proposal["mh_proposal"];
    
    if (Rf_isNewList(proposal_info))
    {
      if ( proposal_info.containsElementNamed("type") && proposal_info.containsElementNamed("variables") )
      {
        std::string type = proposal_info["type"];
        if (type=="norm_rw")
        {
          List info = get_single_variable_one_double_parameter_info(model_parameters,
                                                                    proposal_info,
                                                                    type);
          
          std::string variable = info[0];
          double sd = info[1];
          
          proposal = new GaussianRandomWalkProposalKernel(variable,
                                                          sd);
          
        }
        else if (type=="mvnorm_rw")
        {
          List info = get_single_variable_matrix_parameter_info(model_parameters,
                                                                proposal_info,
                                                                type);
          
          std::string variable = info[0];
          arma::mat cov = info[1];
          
          proposal = new GaussianRandomWalkProposalKernel(variable,
                                                          cov);
          
        }
        if (type=="unif_rw")
        {
          List info = get_single_variable_one_double_parameter_info(model_parameters,
                                                                    proposal_info,
                                                                    type);
          
          std::string variable = info[0];
          double sd = info[1];
          
          proposal = new UniformRandomWalkProposalKernel(variable,
                                                         sd);
          
        }
        if (type=="mvunif_rw")
        {
          List info = get_single_variable_vector_parameter_info(model_parameters,
                                                                proposal_info,
                                                                type);
          
          std::string variable = info[0];
          arma::colvec halfwidth = info[1];
          
          proposal = new UniformRandomWalkProposalKernel(variable,
                                                         halfwidth);
          
        }
        else if (type=="langevin")
        {
          List info = get_single_variable_matrix_parameter_info(model_parameters,
                                                                proposal_info,
                                                                type);
          
          std::string variable = info[0];
          arma::mat cov = info[1];
          
          proposal = new LangevinProposalKernel(variable,
                                                cov,
                                                new DirectGradientEstimator());
          
        }
        else if (type=="hmc")
        {
          List info = get_single_variable_matrix_parameter_info(model_parameters,
                                                                proposal_info,
                                                                type);
          
          std::string variable = info[0];
          arma::mat cov = info[1];
          
          proposal = new HMCProposalKernel(variable,
                                           cov,
                                           new DirectGradientEstimator());
          
        }
        else if (type=="barker_dynamics")
        {
          List info = get_single_variable_matrix_parameter_info(model_parameters,
                                                                proposal_info,
                                                                type);
          
          std::string variable = info[0];
          arma::mat cov = info[1];
          
          proposal = new BarkerDynamicsProposalKernel(variable,
                                                      cov,
                                                      new DirectGradientEstimator());
          
        }
        else if (type=="mirror")
        {
          List info = get_single_variable_vector_and_matrix_parameter_info(model_parameters,
                                                                           proposal_info,
                                                                           type);
          
          std::string variable = info[0];
          arma::colvec mean = info[1];
          arma::mat cov = info[2];
          
          proposal = new MirrorProposalKernel(variable,
                                              mean,
                                              cov);
          
        }
        else
        {
          Rcout << "Proposal type " << type;
          stop("Proposal type unknown");
        }
      }
      else
      {
        stop("Missing information for proposal in model file.");
      }
      
    }
    else
    {
      stop("Error in proposal section of model file.");
    }
  }
  else if ( current_proposal.containsElementNamed("evaluate_log_mh_proposal") && current_proposal.containsElementNamed("simulate_mh_proposal") && current_proposal.containsElementNamed("type") )
  {
    SEXP evaluate_log_mh_proposal_SEXP = current_proposal["evaluate_log_mh_proposal"];
    SEXP simulate_mh_proposal_SEXP = current_proposal["simulate_mh_proposal"];
    size_t type = current_proposal["type"];
    if (type==1)
    {
      proposal = new CustomNoParamsProposalKernel(load_simulate_no_params_mcmc_proposal(simulate_mh_proposal_SEXP),
                                                  load_evaluate_log_no_params_mcmc_proposal(evaluate_log_mh_proposal_SEXP));
    }
    else if (type==2)
    {
      proposal = new CustomProposalKernel(load_simulate_mcmc_proposal(simulate_mh_proposal_SEXP),
                                          load_evaluate_log_mcmc_proposal(evaluate_log_mh_proposal_SEXP));
    }
    else if (type==3)
    {
      proposal = new CustomGuidedNoParamsProposalKernel(load_simulate_guided_no_params_mcmc_proposal(simulate_mh_proposal_SEXP),
                                                        load_evaluate_log_guided_no_params_mcmc_proposal(evaluate_log_mh_proposal_SEXP),data);
    }
    else if (type==4)
    {
      proposal = new CustomGuidedProposalKernel(load_simulate_guided_mcmc_proposal(simulate_mh_proposal_SEXP),
                                                load_evaluate_log_guided_mcmc_proposal(evaluate_log_mh_proposal_SEXP),data);
    }
    else
    {
      stop("Functions read in for proposal are invalid.");
    }
  }
  else
  {
    stop("Invalid mh proposal. Maybe you specified a method to simulate from the proposal, but no method to evaluate it?");
  }
  
  
  if (proposal==NULL)
  {
    stop("No suitable proposal found in model file.");
  }
  else
  {
    return proposal;
  }
}

ProposalKernel* get_unadjusted_proposal(const List &current_proposal,
                                const List &model_parameters,
                                const Data* data)
{
  ProposalKernel* proposal;
  if ( current_proposal.containsElementNamed("unadjusted_proposal") )
  {
    List proposal_info = current_proposal["unadjusted_proposal"];
    
    if (Rf_isNewList(proposal_info))
    {
      if ( proposal_info.containsElementNamed("type") && proposal_info.containsElementNamed("variables") )
      {
        std::string type = proposal_info["type"];
        if (type=="norm_rw")
        {
          List info = get_single_variable_one_double_parameter_info(model_parameters,
                                                                    proposal_info,
                                                                    type);
          
          std::string variable = info[0];
          double sd = info[1];
          
          proposal = new GaussianRandomWalkProposalKernel(variable,
                                                          sd);
          
        }
        else if (type=="mvnorm_rw")
        {
          List info = get_single_variable_matrix_parameter_info(model_parameters,
                                                                proposal_info,
                                                                type);
          
          std::string variable = info[0];
          arma::mat cov = info[1];
          
          proposal = new GaussianRandomWalkProposalKernel(variable,
                                                          cov);
          
        }
        if (type=="unif_rw")
        {
          List info = get_single_variable_one_double_parameter_info(model_parameters,
                                                                    proposal_info,
                                                                    type);
          
          std::string variable = info[0];
          double sd = info[1];
          
          proposal = new UniformRandomWalkProposalKernel(variable,
                                                         sd);
          
        }
        if (type=="mvunif_rw")
        {
          List info = get_single_variable_vector_parameter_info(model_parameters,
                                                                proposal_info,
                                                                type);
          
          std::string variable = info[0];
          arma::colvec halfwidth = info[1];
          
          proposal = new UniformRandomWalkProposalKernel(variable,
                                                         halfwidth);
          
        }
        else if (type=="langevin")
        {
          List info = get_single_variable_matrix_parameter_info(model_parameters,
                                                                proposal_info,
                                                                type);
          
          std::string variable = info[0];
          arma::mat cov = info[1];
          
          proposal = new LangevinProposalKernel(variable,
                                                cov,
                                                new DirectGradientEstimator());
          
        }
        else if (type=="hmc")
        {
          List info = get_single_variable_matrix_parameter_info(model_parameters,
                                                                proposal_info,
                                                                type);
          
          std::string variable = info[0];
          arma::mat cov = info[1];
          
          proposal = new HMCProposalKernel(variable,
                                           cov,
                                           new DirectGradientEstimator());
          
        }
        else if (type=="barker_dynamics")
        {
          List info = get_single_variable_matrix_parameter_info(model_parameters,
                                                                proposal_info,
                                                                type);
          
          std::string variable = info[0];
          arma::mat cov = info[1];
          
          proposal = new BarkerDynamicsProposalKernel(variable,
                                                      cov,
                                                      new DirectGradientEstimator());
          
        }
        else if (type=="mirror")
        {
          List info = get_single_variable_vector_and_matrix_parameter_info(model_parameters,
                                                                           proposal_info,
                                                                           type);
          
          std::string variable = info[0];
          arma::colvec mean = info[1];
          arma::mat cov = info[2];
          
          proposal = new MirrorProposalKernel(variable,
                                              mean,
                                              cov);
          
        }
        else
        {
          Rcout << "Proposal type " << type;
          stop("Proposal type unknown");
        }
      }
      else
      {
        stop("Missing information for proposal in model file.");
      }
      
    }
    else
    {
      stop("Error in proposal section of model file.");
    }
  }
  else if ( current_proposal.containsElementNamed("simulate_unadjusted_proposal") && current_proposal.containsElementNamed("type") )
  {
    SEXP simulate_unadjusted_proposal_SEXP = current_proposal["simulate_unadjusted_proposal"];
    size_t type = current_proposal["type"];
    if (type==1)
    {
      proposal = new CustomNoParamsProposalKernel(load_simulate_no_params_mcmc_proposal(simulate_unadjusted_proposal_SEXP));
    }
    else if (type==2)
    {
      proposal = new CustomProposalKernel(load_simulate_mcmc_proposal(simulate_unadjusted_proposal_SEXP));
    }
    else if (type==3)
    {
      proposal = new CustomGuidedNoParamsProposalKernel(load_simulate_guided_no_params_mcmc_proposal(simulate_unadjusted_proposal_SEXP),
                                                        data);
    }
    else if (type==4)
    {
      proposal = new CustomGuidedProposalKernel(load_simulate_guided_mcmc_proposal(simulate_unadjusted_proposal_SEXP),
                                                data);
    }
    else
    {
      stop("Functions read in for proposal are invalid.");
    }
  }
  else
  {
    stop("Invalid mh proposal. Maybe you specified a method to simulate from the proposal, but no method to evaluate it?");
  }
  
  
  if (proposal==NULL)
  {
    stop("No suitable proposal found in model file.");
  }
  else
  {
    return proposal;
  }
}

IndependentProposalKernel* get_independent_mh_proposal(const List &current_proposal,
                                                       const List &model_parameters,
                                                       const Data* data)
{
  IndependentProposalKernel* proposal;
  if ( current_proposal.containsElementNamed("independent_mh_proposal") )
  {
    List proposal_info = current_proposal["independent_mh_proposal"];
    
    if (Rf_isNewList(proposal_info))
    {
      if ( proposal_info.containsElementNamed("type") && proposal_info.containsElementNamed("variables") )
      {
        std::string type = proposal_info["type"];
        if (type=="norm")
        {
          List info = get_single_variable_two_double_parameter_info(model_parameters,
                                                             proposal_info,
                                                             type);
          
          std::string variable = info[0];
          double mean = info[1];
          double sd = info[2];
          
          proposal = new GaussianIndependentProposalKernel(variable,
                                                           mean,
                                                           sd);
          
        }
        else if (type=="mvnorm")
        {
          List info = get_single_variable_vector_and_matrix_parameter_info(model_parameters,
                                                                           proposal_info,
                                                                           type);
          
          std::string variable = info[0];
          arma::colvec mean = info[1];
          arma::mat cov = info[2];
          
          proposal = new GaussianIndependentProposalKernel(variable,
                                                           mean,
                                                           cov);
          
        }
        else if (type=="lnorm")
        {
          List info = get_single_variable_two_double_parameter_info(model_parameters,
                                                                    proposal_info,
                                                                    type);
          
          std::string variable = info[0];
          double mean = info[1];
          double sd = info[2];
          
          proposal = new LogGaussianIndependentProposalKernel(variable,
                                                              mean,
                                                              sd);
          
        }
        else if (type=="mvlnorm")
        {
          List info = get_single_variable_vector_and_matrix_parameter_info(model_parameters,
                                                                           proposal_info,
                                                                           type);
          
          std::string variable = info[0];
          arma::colvec mean = info[1];
          arma::mat cov = info[2];
          
          proposal = new LogGaussianIndependentProposalKernel(variable,
                                                              mean,
                                                              cov);
          
        }
        else if (type=="gamma")
        {
          List info = get_single_variable_two_double_parameter_info(model_parameters,
                                                             proposal_info,
                                                             type);
          
          std::string variable = info[0];
          double shape = info[1];
          double rate = info[2];
          
          proposal = new GammaIndependentProposalKernel(variable,
                                                        shape,
                                                        rate);
          
        }
        else
        {
          Rcout << "Proposal type " << type;
          stop("Proposal type unknown");
        }
      }
      else
      {
        stop("Missing information for proposal in model file.");
      }
      
    }
    else
    {
      stop("Error in proposal section of model file.");
    }
  }
  else if ( current_proposal.containsElementNamed("evaluate_log_independent_mh_proposal") && current_proposal.containsElementNamed("simulate_independent_mh_proposal") && current_proposal.containsElementNamed("type") )
  {
    SEXP evaluate_log_independent_mh_proposal_SEXP = current_proposal["evaluate_log_independent_mh_proposal"];
    SEXP simulate_independent_mh_proposal_SEXP = current_proposal["simulate_independent_mh_proposal"];
    size_t type = current_proposal["type"];
    if (type==1)
    {
      proposal = new CustomDistributionProposalKernel(load_simulate_distribution(simulate_independent_mh_proposal_SEXP),
                                                      load_evaluate_log_distribution(evaluate_log_independent_mh_proposal_SEXP));
    }
    else if (type==2)
    {
      proposal = new CustomIndependentProposalKernel(load_simulate_independent_proposal(simulate_independent_mh_proposal_SEXP),
                                                     load_evaluate_log_independent_proposal(evaluate_log_independent_mh_proposal_SEXP));
    }
    else if (type==3)
    {
      proposal = new CustomGuidedDistributionProposalKernel(load_simulate_guided_distribution(simulate_independent_mh_proposal_SEXP),
                                                            load_evaluate_log_guided_distribution(evaluate_log_independent_mh_proposal_SEXP),data);
    }
    else if (type==4)
    {
      proposal = new CustomGuidedIndependentProposalKernel(load_simulate_guided_independent_proposal(simulate_independent_mh_proposal_SEXP),
                                                           load_evaluate_log_guided_independent_proposal(evaluate_log_independent_mh_proposal_SEXP),data);
    }
    else
    {
      stop("Functions read in for proposal are invalid.");
    }
  }
  else
  {
    stop("Invalid independent_mh proposal. Maybe you specified a method to simulate from the proposal, but no method to evaluate it?");
  }
  
  if (proposal==NULL)
  {
    stop("No suitable prior to use as proposal.");
  }
  else
  {
    return proposal;
  }
}

SymmetricProposalKernel* get_m_proposal(const List &current_proposal,
                                        const List &model_parameters,
                                        const Data* data)
{
  SymmetricProposalKernel* proposal;
  if ( current_proposal.containsElementNamed("m_proposal") )
  {
    List proposal_info = current_proposal["m_proposal"];
    
    if (Rf_isNewList(proposal_info))
    {
      if ( proposal_info.containsElementNamed("type") && proposal_info.containsElementNamed("variables") )
      {
        std::string type = proposal_info["type"];
        if (type=="norm_rw")
        {
          List info = get_single_variable_one_double_parameter_info(model_parameters,
                                                                    proposal_info,
                                                                    type);
          
          std::string variable = info[0];
          double sd = info[1];
          
          proposal = new GaussianRandomWalkProposalKernel(variable,
                                                          sd);
          
        }
        else if (type=="mvnorm_rw")
        {
          List info = get_single_variable_matrix_parameter_info(model_parameters,
                                                                proposal_info,
                                                                type);
          
          std::string variable = info[0];
          arma::mat cov = info[1];
          
          proposal = new GaussianRandomWalkProposalKernel(variable,
                                                          cov);
          
        }
        if (type=="unif_rw")
        {
          List info = get_single_variable_one_double_parameter_info(model_parameters,
                                                                    proposal_info,
                                                                    type);
          
          std::string variable = info[0];
          double sd = info[1];
          
          proposal = new UniformRandomWalkProposalKernel(variable,
                                                         sd);
          
        }
        if (type=="mvunif_rw")
        {
          List info = get_single_variable_vector_parameter_info(model_parameters,
                                                                proposal_info,
                                                                type);
          
          std::string variable = info[0];
          arma::colvec halfwidth = info[1];
          
          proposal = new UniformRandomWalkProposalKernel(variable,
                                                         halfwidth);
          
        }
        else if (type=="mirror")
        {
          List info = get_single_variable_vector_and_matrix_parameter_info(model_parameters,
                                                                           proposal_info,
                                                                           type);
          
          std::string variable = info[0];
          arma::colvec mean = info[1];
          arma::mat cov = info[2];
          
          proposal = new MirrorProposalKernel(variable,
                                              mean,
                                              cov);
          
        }
        else
        {
          Rcout << "Proposal type " << type;
          stop("Proposal type unknown");
        }
      }
      else
      {
        stop("Missing information for proposal in model file.");
      }
      
    }
    else
    {
      stop("Error in proposal section of model file.");
    }
  }
  else if ( current_proposal.containsElementNamed("simulate_m_proposal") && current_proposal.containsElementNamed("type") )
  {
    SEXP simulate_m_proposal_SEXP = current_proposal["simulate_m_proposal"];
    size_t type = current_proposal["type"];
    if (type==1)
    {
      proposal = new CustomNoParamsSymmetricProposalKernel(load_simulate_no_params_mcmc_proposal(simulate_m_proposal_SEXP));
    }
    else if (type==2)
    {
      proposal = new CustomSymmetricProposalKernel(load_simulate_mcmc_proposal(simulate_m_proposal_SEXP));
    }
    else if (type==3)
    {
      proposal = new CustomGuidedNoParamsSymmetricProposalKernel(load_simulate_guided_no_params_mcmc_proposal(simulate_m_proposal_SEXP),data);
    }
    else if (type==4)
    {
      proposal = new CustomGuidedSymmetricProposalKernel(load_simulate_guided_mcmc_proposal(simulate_m_proposal_SEXP),data);
    }
    else
    {
      stop("Functions read in for proposal are invalid.");
    }
  }
  else
  {
    stop("Invalid m proposal. Maybe you specified a method to simulate from the proposal, but no method to evaluate it?");
  }
  
  if (proposal==NULL)
  {
    stop("No suitable proposal found in model file.");
  }
  else
  {
    return proposal;
  }
}

MCMCTermination* make_mcmc_termination(const List &model_parameters,
                                       const List &mcmc_termination_method)
{
  MCMCTermination* termination_method;
  
  if (mcmc_termination_method.containsElementNamed("method") && mcmc_termination_method.containsElementNamed("values"))
  {
    std::string method = mcmc_termination_method["method"];
    if (method=="iterations")
    {
      if (Rf_isNewList(mcmc_termination_method["values"]))
      {
        List values = mcmc_termination_method["values"];
        
        if (values.length()==1)
        {
          int number_of_iterations = extract_int_parameter(values,
                                                           model_parameters,
                                                           0);
          termination_method = new IterationsMCMCTermination(number_of_iterations);
        }
        else
        {
          stop("MCMC termination using iterations requires you to specify only the number of iterations.");
        }
      }
      else
      {
        stop("MCMC termination using iterations requires you to specify only the number of iterations.");
      }
    }
    else if (method=="se")
    {
      if (Rf_isNewList(mcmc_termination_method["values"]))
      {
        List values = mcmc_termination_method["values"];
        
        if (values.length()==1)
        {
          double threshold = extract_double_parameter(values,
                                                      model_parameters,
                                                      0);
          int number_of_iterations = arma::datum::inf;
          
          termination_method = new SEMCMCTermination(threshold,
                                                     number_of_iterations);
        }
        else if (values.length()==2)
        {
          double threshold = extract_double_parameter(values,
                                                      model_parameters,
                                                      0);
          int number_of_iterations = extract_int_parameter(values,
                                                           model_parameters,
                                                           1);
          
          termination_method = new SEMCMCTermination(threshold,
                                                     number_of_iterations);
        }
        else
        {
          stop("MCMC termination using standard error requires you to specify a threshold.");
        }
      }
      else
      {
        stop("MCMC termination using standard error requires you to specify a threshold.");
      }
    }
    else
    {
      stop("Invalid method for MCMC termination.");
    }
  }
  else
  {
    stop("No valid method found for MCMC termination.");
  }
  
  return termination_method;
}

MCMC* make_mcmc(const List &model,
                const List &model_parameters,
                const Data* data,
                const List &mcmc_termination_method,
                const List &mcmc_weights_method,
                Index* full_index)
{
  MCMCTermination* termination_method = make_mcmc_termination(model_parameters,mcmc_termination_method);
  
  MCMC* mcmc = NULL;
  
  if (model.containsElementNamed("order_of_mcmc"))
  {
    NumericVector order_of_mcmc = model["order_of_mcmc"];
    std::vector<MCMC*> moves;
    moves.reserve(order_of_mcmc.size());
    
    size_t simulate_mh_proposal_index = 0;
    size_t simulate_unadjusted_proposal_index = 0;
    size_t simulate_independent_mh_proposal_index = 0;
    size_t simulate_m_proposal_index = 0;
    
    for (auto i=order_of_mcmc.begin();
         i!=order_of_mcmc.end();
         ++i)
    {
      
      if ((*i)==1)
      {
        ProposalKernel* proposal;
        
        if (model.containsElementNamed("mh_proposal"))
        {
          List proposal_infos = model["mh_proposal"];
          if (Rf_isNewList(proposal_infos[simulate_mh_proposal_index]))
          {
            proposal = get_mh_proposal(proposal_infos[simulate_mh_proposal_index],
                                       model_parameters,
                                       data);
          }
          else
          {
            stop("Error in mcmc proposals part of model file.");
          }
          
          if (order_of_mcmc.size()>1)
          {
            mcmc = new MetropolisHastingsMCMC(1,
                                              proposal);
          }
          else
          {
            mcmc = new MetropolisHastingsMCMC(termination_method,
                                              proposal);
          }
          
          if (proposal_infos.containsElementNamed("mh_factor_index"))
          {
            SEXP index_SEXP = proposal_infos["mh_factor_index"];
            NumericVector index = load_numeric_vector(index_SEXP);
            std::vector<size_t> indices_in;
            for (size_t k=0; k<index.length(); ++k)
            {
              indices_in.push_back(size_t(index[k]));
            }
            mcmc->set_index_if_null(new VectorSingleIndex(indices_in));
          }
          
          simulate_mh_proposal_index = simulate_mh_proposal_index + 1;
        }
        else
        {
          stop("No mh_proposal found.");
        }
        
      }
      else if ((*i)==2)
      {
        IndependentProposalKernel* proposal;
        
        if (model.containsElementNamed("independent_mh_proposal"))
        {
          List proposal_infos = model["independent_mh_proposal"];
          
          if (Rf_isNewList(proposal_infos[simulate_independent_mh_proposal_index]))
          {
            proposal = get_independent_mh_proposal(proposal_infos[simulate_independent_mh_proposal_index],
                                                   model_parameters,
                                                   data);
          }
          else
          {
            stop("Error in mcmc proposals part of model file.");
          }
          
          if (order_of_mcmc.size()>1)
          {
            mcmc = new MetropolisHastingsMCMC(1,
                                              proposal);
          }
          else
          {
            mcmc = new MetropolisHastingsMCMC(termination_method,
                                              proposal);
          }
          
          if (proposal_infos.containsElementNamed("independent_mh_factor_index"))
          {
            SEXP index_SEXP = proposal_infos["independent_mh_factor_index"];
            NumericVector index = load_numeric_vector(index_SEXP);
            std::vector<size_t> indices_in;
            for (size_t k=0; k<index.length(); ++k)
            {
              indices_in.push_back(size_t(index[k]));
            }
            mcmc->set_index_if_null(new VectorSingleIndex(indices_in));
          }
          
          simulate_independent_mh_proposal_index = simulate_independent_mh_proposal_index + 1;
        }
        else
        {
          stop("No independent_mh_proposal found.");
        }
        
      }
      else if ((*i)==3)
      {
        SymmetricProposalKernel* proposal;
        
        if (model.containsElementNamed("m_proposal"))
        {
          List proposal_infos = model["m_proposal"];
          
          if (Rf_isNewList(proposal_infos[simulate_m_proposal_index]))
          {
            proposal = get_m_proposal(proposal_infos[simulate_m_proposal_index],
                                      model_parameters,
                                      data);
          }
          else
          {
            stop("Error in mcmc proposals part of model file.");
          }
          
          if (order_of_mcmc.size()>1)
          {
            mcmc = new MetropolisMCMC(1,
                                      proposal);
          }
          else
          {
            mcmc = new MetropolisMCMC(termination_method,
                                      proposal);
          }
          
          if (proposal_infos.containsElementNamed("m_factor_index"))
          {
            SEXP index_SEXP = proposal_infos["m_factor_index"];
            NumericVector index = load_numeric_vector(index_SEXP);
            std::vector<size_t> indices_in;
            for (size_t k=0; k<index.length(); ++k)
            {
              indices_in.push_back(size_t(index[k]));
            }
            mcmc->set_index_if_null(new VectorSingleIndex(indices_in));
          }
          
          simulate_m_proposal_index = simulate_m_proposal_index + 1;
        }
        else
        {
          stop("No m_proposal found.");
        }
        
      }
      else if ((*i)==4)
      {
        ProposalKernel* proposal;
        
        if (model.containsElementNamed("unadjusted_proposal"))
        {
          List proposal_infos = model["unadjusted_proposal"];
          if (Rf_isNewList(proposal_infos[simulate_unadjusted_proposal_index]))
          {
            proposal = get_unadjusted_proposal(proposal_infos[simulate_unadjusted_proposal_index],
                                               model_parameters,
                                               data);
          }
          else
          {
            stop("Error in mcmc proposals part of model file.");
          }
          
          if (order_of_mcmc.size()>1)
          {
            mcmc = new UnadjustedMCMC(1,
                                      proposal);
          }
          else
          {
            mcmc = new UnadjustedMCMC(termination_method,
                                      proposal);
          }
          
          if (proposal_infos.containsElementNamed("unadjusted_factor_index"))
          {
            SEXP index_SEXP = proposal_infos["unadjusted_factor_index"];
            NumericVector index = load_numeric_vector(index_SEXP);
            std::vector<size_t> indices_in;
            for (size_t k=0; k<index.length(); ++k)
            {
              indices_in.push_back(size_t(index[k]));
            }
            mcmc->set_index_if_null(new VectorSingleIndex(indices_in));
          }
          
          simulate_unadjusted_proposal_index = simulate_unadjusted_proposal_index + 1;
        }
        else
        {
          stop("No unadjusted_proposal found.");
        }
        
      }
      else
      {
        Rcpp::stop("get_mcmc: invalid type in order_of_mcmc.");
      }
      
      moves.push_back(mcmc);
      
    }
    
    if (mcmc_weights_method.containsElementNamed("func"))
    {
      NumericVector mcmc_weights;

      SEXP mcmc_weights_SEXP = mcmc_weights_method["func"];
      mcmc_weights = load_numeric_vector(mcmc_weights_SEXP);
      
      if (moves.size()==1)
      {
        mcmc->set_index_if_null(full_index);
        return mcmc;
      }
      else
      {
        MCMC* new_mcmc = new StochasticScanMCMC(termination_method,
                                                moves,
                                                mcmc_weights);
        new_mcmc->set_index_if_null(full_index);
        return new_mcmc;
      }
      
    }
    else
    {
      if (moves.size()==1)
      {
        mcmc->set_index_if_null(full_index);
        return mcmc;
      }
      else
      {
        MCMC* new_mcmc = new DeterministicScanMCMC(termination_method,
                                                   moves);
        new_mcmc->set_index_if_null(full_index);
        return new_mcmc;
      }
    }
  }
  else
  {
    Rcpp::stop("get_mcmc: order_of_mcmc not specified in model.");
  }
}

std::vector<Parameters> make_initial_points(const List &initial_values)
{
  std::vector<Parameters> output;
  output.reserve(initial_values.size());
  
  for (size_t i=0; i<initial_values.size(); ++i)
  {
    Parameters current_parameters;
    if (Rf_isNewList(initial_values[i]))
    {
      List list_params = initial_values[i];
      if (list_params.size()>0)
      {
        CharacterVector names = list_params.names();
        
        for (size_t i=0; i<names.size(); ++i)
        {
          if (isDoubleInList(list_params,Rcpp::as<std::string>(names[i])))
          {
            double param = list_params[Rcpp::as<std::string>(names[i])];
            current_parameters[Rcpp::as<std::string>(names[i])] = param;
          }
          else if (isNumericVectorInList(list_params,Rcpp::as<std::string>(names[i])))
          {
            NumericVector param = list_params[Rcpp::as<std::string>(names[i])];
            current_parameters[Rcpp::as<std::string>(names[i])] = as<arma::colvec>(param);
          }
          else if (isNumericMatrixInList(list_params,Rcpp::as<std::string>(names[i])))
          {
            arma::mat param = list_params[Rcpp::as<std::string>(names[i])];
            current_parameters[Rcpp::as<std::string>(names[i])] = param;
          }
          else
          {
            stop("Initial value must be a float/double, vector or matrix.");
          }
        }
      }
      else
      {
        Rcpp::stop("List of initial values has empty entries.");
      }
      
      output.push_back(current_parameters);
    }
    else
    {
      Rcpp::stop("List of initial values must be a list (over chains) of lists, each of which contain the initial parameters.");
    }
  }

  return output;
}

// [[Rcpp::export]]
size_t ilike_rdtsc()
{
  return rdtsc();
}
  
Parameters make_algorithm_parameters(const List &algorithm_parameter_list)
{
  /*
  Parameters output;
  if (algorithm_parameter_list.size()>0)
  {
    CharacterVector names = algorithm_parameter_list.names();
    
    for (size_t i=0; i<names.size(); ++i)
    {
      arma::mat param = algorithm_parameter_list[Rcpp::as<std::string>(names[i])];
      output[Rcpp::as<std::string>(names[i])] = param;
    }
  }
  */
  
  return list_to_parameters(algorithm_parameter_list);
}

SMCCriterion* get_resampling_method(const List &model,
                                    const List &model_parameters,
                                    const List &adaptive_resampling_method,
                                    size_t number_of_particles)
{
  SMCCriterion* smc_method;
  
  if (adaptive_resampling_method.containsElementNamed("method") && adaptive_resampling_method.containsElementNamed("values"))
  {
    std::string method = adaptive_resampling_method["method"];
    if (method=="ess")
    {
      if (Rf_isNewList(adaptive_resampling_method["values"]))
      {
        List values = adaptive_resampling_method["values"];
        
        if (values.length()==1)
        {
          double proportion = extract_double_parameter(values,
                                                       model_parameters,
                                                       0);
          smc_method = new ESSSMCCriterion(proportion*double(number_of_particles));
        }
        else
        {
          stop("Adaptive resampling using ESS requires specification of the proportion of the total number of particles.");
        }
      }
      else
      {
        stop("Adaptive resampling using ESS requires specification of the proportion of the total number of particles.");
      }
    }
    else
    {
      stop("No valid method found for adaptive resampling.");
    }
  }
  else
  {
    stop("No valid method found for adaptive resampling (need method and values).");
  }
  
  return smc_method;
}

SMCCriterion* get_adaptive_target_method(const List &model,
                                         const List &model_parameters,
                                         const List &adaptive_target_method,
                                         size_t number_of_particles)
{
  SMCCriterion* smc_method;
  
  if (adaptive_target_method.containsElementNamed("method"))
  {
    std::string method = adaptive_target_method["method"];
    if (method=="ess")
    {
      if (adaptive_target_method.containsElementNamed("values"))
      {
        if (Rf_isNewList(adaptive_target_method["values"]))
        {
          List values = adaptive_target_method["values"];
          
          if (values.length()==1)
          {
            double proportion = extract_double_parameter(values,
                                                         model_parameters,
                                                         0);
            smc_method = new ESSSMCCriterion(proportion*double(number_of_particles));
          }
          else
          {
            stop("Adaptive target using ESS requires specification of the proportion of the total number of particles.");
          }
        }
        else
        {
          stop("No valid method found for adaptive target.");
        }
      }
      else
      {
        stop("Adaptive target using ESS requires specification of the proportion of the total number of particles (values need to be specified).");
      }
    }
    else if (method=="cess")
    {
      if (adaptive_target_method.containsElementNamed("values"))
      {
        if (Rf_isNewList(adaptive_target_method["values"]))
        {
          List values = adaptive_target_method["values"];
          
          if (values.length()==1)
          {
            double proportion = extract_double_parameter(values,
                                                         model_parameters,
                                                         0);
            //Rcout << proportion << std::endl;
            //Rcout << proportion*double(number_of_particles) << std::endl;
            smc_method = new CESSSMCCriterion(proportion*double(number_of_particles));
          }
          else
          {
            stop("Adaptive target using CESS requires specification of the proportion of the total number of particles.");
          }
        }
        else
        {
          stop("No valid method found for adaptive target.");
        }
      }
      else
      {
        stop("Adaptive target using ESS requires specification the proportion of the total number of particles (values need to be specified).");
      }
    }
    else if (method=="positive")
    {
      smc_method = new PositiveSMCCriterion();
    }
    else
    {
      stop("No valid method found for adaptive target.");
    }
  }
  else
  {
    stop("No valid method found for adaptive target (need method).");
  }
  
  return smc_method;
}

SMCTermination* get_smc_termination_method(const List &model,
                                           const List &model_parameters,
                                           const List &smc_termination_method)
{
  SMCTermination* smc_termination;
  
  if (smc_termination_method.length()==0)
  {
    smc_termination = NULL;
    return smc_termination;
  }
  
  if (smc_termination_method.containsElementNamed("method"))
  {
    std::string method = smc_termination_method["method"];
    if (method=="stable")
    {
      if (smc_termination_method.containsElementNamed("values"))
      {
        if (Rf_isNewList(smc_termination_method["values"]))
        {
          List values = smc_termination_method["values"];
          
          if (values.length()==2)
          {
            double threshold = extract_double_parameter(values,
                                                        model_parameters,
                                                        0);
            double number_in_a_row = extract_int_parameter(values,
                                                           model_parameters,
                                                           1);
            smc_termination = new StableSMCTermination(threshold,
                                                       number_in_a_row);
          }
          else
          {
            stop("SMC termination via stability requires the specificaton of a threshold and the number of iterations in a row that this threshold must be exceeded.");
          }
        }
        else
        {
          stop("No valid method found for SMC termination.");
        }
      }
      else
      {
        stop("SMC termination via stability requires specification of the desired ESS (values need to be specified).");
      }
      
    }
    else if (method=="always")
    {
      smc_termination = new AlwaysSMCTermination();
    }
    else
    {
      stop("No valid method found for SMC termination.");
    }
  }
  else
  {
    stop("No valid method found for SMC termination (method needed).");
  }
  
  return smc_termination;
}

EnsembleShifter* get_enk_shifter_method(const List &model,
                                    const List &model_parameters,
                                    const List &enk_shifter_method)
{
  EnsembleShifter* enk_shifter;
  
  if (enk_shifter_method.length()==0)
  {
    enk_shifter = NULL;
    return enk_shifter;
  }
  
  if (enk_shifter_method.containsElementNamed("method"))
  {
    std::string method = enk_shifter_method["method"];
    if (method=="stochastic")
    {
      enk_shifter = new StochasticEnsembleShifter();
    }
    else
    {
      stop("No valid method found for EnK shifter.");
    }
  }
  else
  {
    stop("No valid method found for EnK shifter (method needed).");
  }
  
  return enk_shifter;
}

List get_smc_sequencer_info(const List &model,
                            const List &model_parameters,
                            const List &smc_sequencer_method)
{
  if (smc_sequencer_method.containsElementNamed("types") && smc_sequencer_method.containsElementNamed("variables") && smc_sequencer_method.containsElementNamed("schedules"))
  {
    int number_of_bisections;
    if (smc_sequencer_method.containsElementNamed("bisections"))
    {
      List values = smc_sequencer_method["bisections"];
      number_of_bisections = extract_int_parameter(values,
                                                   model_parameters,
                                                   0);
    }
    else
    {
      Rcout << "Bisections not set in adaptive target method, using 100." << std::endl;
      number_of_bisections = 100;
    }
    
    StringVector types;
    if (isStringVectorInList(smc_sequencer_method, "types"))
    {
      types = smc_sequencer_method["types"];
    }
    
    StringVector variables;
    if (isStringVectorInList(smc_sequencer_method, "variables"))
    {
      variables = smc_sequencer_method["variables"];
    }
    
    List schedules;
    if (Rf_isNewList(smc_sequencer_method["schedules"]))
    {
      schedules = smc_sequencer_method["schedules"];
    }
    else if (isNumericVectorInList(smc_sequencer_method,"schedules"))
    {
      schedules = List::create(smc_sequencer_method["schedules"]);
    }
    else
    {
      stop("Invalid specification of schedules for SMC: must be either a numeric vector or a list of numeric vectors.");
    }
    
    if ( (types.size()==variables.size()) && (types.size()==schedules.size()) )
    {
      std::vector<std::string> types_for_output;
      types_for_output.reserve(types.size());
      for (size_t i=0; i<types.size(); ++i)
      {
        types_for_output.push_back(String(types[i]).get_cstring());
      }
      
      std::vector<std::string> variables_for_output;
      variables_for_output.reserve(variables.size());
      for (size_t i=0; i<variables.size(); ++i)
      {
        variables_for_output.push_back(String(variables[i]).get_cstring());
      }
      
      std::vector<std::vector<double>> schedules_for_output;
      schedules_for_output.reserve(schedules.size());
      for (size_t i=0; i<schedules.size(); ++i)
      {
        std::vector<double> schedule_for_output;
        
        /*
        arma::colvec schedule = extract_vector_parameter(schedules,
        model_parameters,
        i);
        schedule_for_output = arma::conv_to<std::vector<double>>::from(schedule);
        */
        
        if (isNumericVectorInList(schedules,i+1))
        {
          NumericVector schedule = schedules[i];
          schedule_for_output = as<std::vector<double>>(schedule);
        }
        else
        {
          stop("Schedule invalid: needs to be a numeric vector.");
        }
        
        schedules_for_output.push_back(schedule_for_output);
      }
      
      if (smc_sequencer_method.containsElementNamed("factors_affected"))
      {
        if (isNumericVectorInList(smc_sequencer_method,"factors_affected"))
        {
          NumericVector factors_affected = smc_sequencer_method["factors_affected"];
          //std::vector<int> factors_affected_for_output = as<std::vector<int>>(factors_affected);
          
          return List::create(types_for_output,variables_for_output,schedules_for_output,number_of_bisections,factors_affected);
        }
        else
        {
          stop("Factors affected invalid: needs to be a numeric vector.");
        }
      }
      else
      {
        return List::create(types_for_output,variables_for_output,schedules_for_output,number_of_bisections);
      }
    }
    else
    {
      stop("Invalid specification of schedules for SMC: types, variables and schedules (list) must all be the same length.");
    }
    
  }
  else
  {
    std::vector<std::string> types_for_output;
    std::vector<std::string> variables_for_output;
    std::vector<std::vector<double>> schedules_for_output;
    size_t number_of_bisections = 100;
    return List::create(types_for_output,variables_for_output,schedules_for_output,number_of_bisections);
    //stop("No valid method found for SMC sequence: need to specify types, variables and schedules.");
  }
}

/*
std::vector<LikelihoodEstimator*> convent_to_annealed_likelihoods_if_needed(RandomNumberGenerator* rng_in,
                                                                            size_t* seed_in,
                                                                            Data* data_in,
                                                                            const std::vector<std::string> &sequencer_types,
                                                                            const std::vector<std::string> &sequencer_variables,
                                                                            const std::vector<LikelihoodEstimator*> &likelihood_estimators,
                                                                            IndependentProposalKernel* proposal,
                                                                            bool proposal_is_evaluated_in)
{
  bool any_annealing = false;
  std::string annealing_variable;
  for (size_t i=0; i<sequencer_types.size(); ++i)
  {
    if ( (sequencer_types[i]=="annealing") || (sequencer_types[i]=="tempering") )
    {
      if (any_annealing==false)
      {
        any_annealing = true;
        annealing_variable = sequencer_variables[i];
        break;
      }
      else
      {
        stop("Two annealing variables specified: not currently supported.");
      }
    }
  }
  
  if (!any_annealing)
  {
    return likelihood_estimators;
  }
  
  std::vector<LikelihoodEstimator*> annealed_likelihood_estimators;
  annealed_likelihood_estimators.reserve(likelihood_estimators.size());
  
  PowerFunctionPtr power = annealing_power;
  
  for (auto i=likelihood_estimators.begin();
       i!=likelihood_estimators.end();
       ++i)
  {
    annealed_likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(rng_in,
                                                                             seed_in,
                                                                             data_in,
                                                                             *i,
                                                                             power,
                                                                             annealing_variable,
                                                                             false));
  }
  
  if (proposal_is_evaluated_in)
  {
    ExactLikelihoodEstimator* proposal_for_evaluation = new ExactLikelihoodEstimator(rng_in,
                                                                                     seed_in,
                                                                                     data_in,
                                                                                     proposal,
                                                                                     true);
    
    PowerFunctionPtr second_power = annealing_one_minus_power;
    
    annealed_likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(rng_in,
                                                                             seed_in,
                                                                             data_in,
                                                                             proposal_for_evaluation,
                                                                             second_power,
                                                                             annealing_variable,
                                                                             false));
    
  }

  return annealed_likelihood_estimators;
}
*/

// [[Rcpp::export]]
double do_importance_sampler(const List &model,
                             const List &parameters,
                             const List &algorithm_parameter_list,
                             size_t number_of_importance_points,
                             bool parallel_in,
                             size_t grain_size_in,
                             const String &results_name_in,
                             size_t seed)
{
  RandomNumberGenerator rng;
  Data the_data = get_data(model);
  std::vector<Data> data_created_in_get_likelihood_estimators;
  
  //std::string results_name = "/Users/richard/Dropbox/code/ilike/experiments/test";
  
  // May need to alter for cases where the likelihood needs to be tuned automatically (e.g. in ABC).
  
  // Check if the prior is the proposal: affects what llhd_estimators we include.
  
  std::vector<LikelihoodEstimator*> likelihood_estimators;
  IndependentProposalKernel* proposal_in;
  bool proposal_is_evaluated_in;
  
  List sequencer_info = get_smc_sequencer_info(model,
                                               parameters,
                                               List());
  std::vector<std::string> sequencer_types = sequencer_info[0];
  std::vector<std::string> sequencer_variables = sequencer_info[1];
  std::vector<std::vector<double>> sequencer_schedules = sequencer_info[2];
  
  std::vector<int> factors_affected_by_smc_sequence;
  if (sequencer_info.length()==5)
  {
    NumericVector factors_affected_by_smc_sequence_temp = sequencer_info[4];
    factors_affected_by_smc_sequence = as<std::vector<int>>(factors_affected_by_smc_sequence_temp);
  }
  
  Index* without_cancelled_index = NULL;
  Index* full_index = NULL;
  bool any_annealing = false;
  
  
  if ( model.containsElementNamed("importance_proposal") )
  {
    
    proposal_in = get_proposal(model,
                               parameters,
                               &the_data);
    
    likelihood_estimators = get_likelihood_estimators(&rng,
                                                      &seed,
                                                      &the_data,
                                                      model,
                                                      parameters,
                                                      true,
                                                      sequencer_types,
                                                      sequencer_variables,
                                                      sequencer_schedules,
                                                      proposal_in,
                                                      without_cancelled_index,
                                                      full_index,
                                                      any_annealing,
                                                      factors_affected_by_smc_sequence,
                                                      data_created_in_get_likelihood_estimators);
    
    proposal_is_evaluated_in = true;
    
  }
  else
  {
    Rcout << "No importance sampling proposal specified; using prior." << std::endl;
    
    proposal_in = get_prior_as_simulate_only_proposal(model,
                                                      parameters);
    
    likelihood_estimators = get_likelihood_estimators(&rng,
                                                      &seed,
                                                      &the_data,
                                                      model,
                                                      parameters,
                                                      false,
                                                      sequencer_types,
                                                      sequencer_variables,
                                                      sequencer_schedules,
                                                      proposal_in,
                                                      without_cancelled_index,
                                                      full_index,
                                                      any_annealing,
                                                      factors_affected_by_smc_sequence,
                                                      data_created_in_get_likelihood_estimators);
    
    proposal_is_evaluated_in = false;
  }
  
  
  Parameters algorithm_parameters = make_algorithm_parameters(algorithm_parameter_list);
  
  ImportanceSampler alg(&rng,
                        &seed,
                        &the_data,
                        algorithm_parameters,
                        number_of_importance_points,
                        "",
                        likelihood_estimators,
                        proposal_in,
                        proposal_is_evaluated_in,
                        true,
                        true,
                        false,
                        parallel_in,
                        grain_size_in,
                        "");
  
  
  
  SMCOutput* output = alg.run();
  
  if (strcmp(results_name_in.get_cstring(),"") != 0)
    output->write(results_name_in.get_cstring());
  
  double log_likelihood = output->log_likelihood;
  
  delete output;
  
  return log_likelihood;
}
  
// [[Rcpp::export]]
void do_mcmc(const List &model,
             const List &parameters,
             const List &algorithm_parameter_list,
             const List &initial_values,
             const List &mcmc_termination_method,
             const List &mcmc_weights_method,
             size_t number_of_chains,
             bool parallel_in,
             size_t grain_size_in,
             const String &results_name_in,
             size_t seed)
{
  // could specify with mhproposalkernels
  // if more than one, need to specify type of sweep
  // for each need to specify if g, m or mh proposal, also which factors involved
  
  // either use initial values list if supplied, else use prior
  
  // ilike:: includes any of the options we have made
  
  RandomNumberGenerator rng;
  Data the_data = get_data(model);
  std::vector<Data> data_created_in_get_likelihood_estimators;
  
  // May need to alter for cases where the likelihood needs to be tuned automatically (e.g. in ABC).
  
  // Check if the prior is the proposal: affects what llhd_estimators we include.
  
  List sequencer_info = get_smc_sequencer_info(model,
                                               parameters,
                                               List());
  std::vector<std::string> sequencer_types = sequencer_info[0];
  std::vector<std::string> sequencer_variables = sequencer_info[1];
  std::vector<std::vector<double>> sequencer_schedules = sequencer_info[2];
  
  std::vector<int> factors_affected_by_smc_sequence;
  if (sequencer_info.length()==5)
  {
    NumericVector factors_affected_by_smc_sequence_temp = sequencer_info[4];
    factors_affected_by_smc_sequence = as<std::vector<int>>(factors_affected_by_smc_sequence_temp);
  }
  
  Index* without_cancelled_index = NULL;
  Index* full_index = NULL;
  bool any_annealing = false;
  
  std::vector<LikelihoodEstimator*> likelihood_estimators;
  likelihood_estimators = get_likelihood_estimators(&rng,
                                                    &seed,
                                                    &the_data,
                                                    model,
                                                    parameters,
                                                    true,
                                                    sequencer_types,
                                                    sequencer_variables,
                                                    sequencer_schedules,
                                                    NULL,
                                                    without_cancelled_index,
                                                    full_index,
                                                    any_annealing,
                                                    factors_affected_by_smc_sequence,
                                                    data_created_in_get_likelihood_estimators);
  
  Parameters algorithm_parameters = make_algorithm_parameters(algorithm_parameter_list);
  
  MCMC* the_mcmc = make_mcmc(model,
                             parameters,
                             &the_data,
                             mcmc_termination_method,
                             mcmc_weights_method,
                             full_index);
  
  SMCMCMCMove* alg;
  
  
  if (initial_values.size()==0)
  {
    
    IndependentProposalKernel* proposal_in = get_prior_as_proposal(model,
                                                                   parameters);
    
    /*
    alg = new SMCMCMCMove(&rng,
                          &seed,
                          &the_data,
                          algorithm_parameters,
                          number_of_chains,
                          2,
                          2,
                          the_mcmc,
                          likelihood_estimators,
                          proposal_in,
                          without_cancelled_index,
                          full_index,
                          true,
                          parallel_in,
                          grain_size_in,
                          "");
    */
    
    alg = new SMCMCMCMove(&rng,
                          &seed,
                          &the_data,
                          algorithm_parameters,
                          number_of_chains,
                          2,
                          2,
                          the_mcmc,
                          likelihood_estimators,
                          proposal_in,
                          without_cancelled_index,
                          full_index,
                          true,
                          parallel_in,
                          grain_size_in,
                          "");
  }
  else
  {
    std::vector<Parameters> initial_points = make_initial_points(initial_values);
    
    
    arma::colvec log_probabilities_of_initial_values(initial_points.size());
    
    alg = new SMCMCMCMove(&rng,
                          &seed,
                          &the_data,
                          algorithm_parameters,
                          2,
                          2,
                          the_mcmc,
                          likelihood_estimators,
                          initial_points,
                          log_probabilities_of_initial_values,
                          without_cancelled_index,
                          full_index,
                          true,
                          parallel_in,
                          grain_size_in,
                          "");
    
  }
  
  
  SMCOutput* output = alg->run();

  
  if (strcmp(results_name_in.get_cstring(),"") != 0)
    output->write(results_name_in.get_cstring());
  
  
  delete output;
  
  delete alg;
  

}


// [[Rcpp::export]]
double do_smc_mcmc_move(const List &model,
                        const List &parameters,
                        const List &algorithm_parameter_list,
                        size_t number_of_particles,
                        const List &mcmc_termination_method,
                        const List &mcmc_weights_method,
                        const List &adaptive_resampling_method,
                        const List &smc_sequencer_method,
                        const List &adaptive_target_method,
                        const List &smc_termination_method,
                        size_t smc_iterations_to_store,
                        bool write_to_file_at_each_iteration,
                        bool parallel_in,
                        size_t grain_size_in,
                        const String &results_name_in,
                        size_t seed)
{
  
  RandomNumberGenerator rng;
  
  Data the_data = get_data(model);
  std::vector<Data> data_created_in_get_likelihood_estimators;
  
  //std::string results_name = "/Users/richard/Dropbox/code/ilike/experiments/test";
  
  // May need to alter for cases where the likelihood needs to be tuned automatically (e.g. in ABC).
  
  // Check if the prior is the proposal: affects what llhd_estimators we include.
  
  std::vector<LikelihoodEstimator*> likelihood_estimators;
  IndependentProposalKernel* proposal_in;
  bool proposal_is_evaluated_in;
  
  List sequencer_info = get_smc_sequencer_info(model,
                                               parameters,
                                               smc_sequencer_method);
  std::vector<std::string> sequencer_types = sequencer_info[0];
  std::vector<std::string> sequencer_variables = sequencer_info[1];
  std::vector<std::vector<double>> sequencer_schedules = sequencer_info[2];
  
  std::vector<int> factors_affected_by_smc_sequence;
  if (sequencer_info.length()==5)
  {
    NumericVector factors_affected_by_smc_sequence_temp = sequencer_info[4];
    factors_affected_by_smc_sequence = as<std::vector<int>>(factors_affected_by_smc_sequence_temp);
  }
  
  Index* without_cancelled_index = NULL;
  Index* full_index = NULL;
  
  bool any_annealing = false;
  
  if ( model.containsElementNamed("importance_proposal") )
  {
    proposal_in = get_proposal(model,
                               parameters,
                               &the_data);
    
    likelihood_estimators = get_likelihood_estimators(&rng,
                                                      &seed,
                                                      &the_data,
                                                      model,
                                                      parameters,
                                                      true,
                                                      sequencer_types,
                                                      sequencer_variables,
                                                      sequencer_schedules,
                                                      proposal_in,
                                                      without_cancelled_index,
                                                      full_index,
                                                      any_annealing,
                                                      factors_affected_by_smc_sequence,
                                                      data_created_in_get_likelihood_estimators);
    
    proposal_is_evaluated_in = true;
    
  }
  else
  {
    Rcout << "No importance sampling proposal specified; using prior." << std::endl;
    
    proposal_in = get_prior_as_simulate_only_proposal(model,
                                                      parameters);
    
    likelihood_estimators = get_likelihood_estimators(&rng,
                                                      &seed,
                                                      &the_data,
                                                      model,
                                                      parameters,
                                                      false,
                                                      sequencer_types,
                                                      sequencer_variables,
                                                      sequencer_schedules,
                                                      proposal_in,
                                                      without_cancelled_index,
                                                      full_index,
                                                      any_annealing,
                                                      factors_affected_by_smc_sequence,
                                                      data_created_in_get_likelihood_estimators);
    
    proposal_is_evaluated_in = false;
    
  }
  
  // note that this will always result in the temperature being adapted in the SMC algorithm (not smcfixed)
  /*
  likelihood_estimators = convent_to_annealed_likelihoods_if_needed(&rng,
                                                                    &seed,
                                                                    &the_data,
                                                                    sequencer_info[0],
                                                                    sequencer_info[1],
                                                                    likelihood_estimators,
                                                                    proposal_in,
                                                                    proposal_is_evaluated_in);
  */
  
  Parameters algorithm_parameters = make_algorithm_parameters(algorithm_parameter_list);
  
  MCMC* the_mcmc = make_mcmc(model,
                             parameters,
                             &the_data,
                             mcmc_termination_method,
                             mcmc_weights_method,
                             full_index);
  
  SMCCriterion* resampling_criterion = get_resampling_method(model,
                                                             parameters,
                                                             adaptive_resampling_method,
                                                             number_of_particles);
  
  SMCCriterion* adaptive_target_criterion = get_adaptive_target_method(model,
                                                                       parameters,
                                                                       adaptive_target_method,
                                                                       number_of_particles);
  
  SMCTermination* smc_termination = get_smc_termination_method(model,
                                                               parameters,
                                                               smc_termination_method);
  
  /*
  if (without_cancelled_index==NULL)
  {
    Rcout << "Without cancelled index is NULL." << std::endl;
  }
  
  if (full_index==NULL)
  {
    Rcout << "Full index is NULL." << std::endl;
  }
  */
  
  if (write_to_file_at_each_iteration)
  {
    SMCMCMCMove* alg = new SMCMCMCMove(&rng,
                                       &seed,
                                       &the_data,
                                       algorithm_parameters,
                                       number_of_particles,
                                       smc_iterations_to_store,
                                       smc_iterations_to_store,
                                       the_mcmc,
                                       resampling_criterion,
                                       adaptive_target_criterion,
                                       sequencer_info[3],
                                       smc_termination,
                                       sequencer_variables,
                                       sequencer_schedules,
                                       likelihood_estimators,
                                       proposal_in,
                                       without_cancelled_index,
                                       full_index,
                                       proposal_is_evaluated_in,
                                       true,
                                       true,
                                       false,
                                       false,
                                       parallel_in,
                                       grain_size_in,
                                       results_name_in.get_cstring());
    
    //Rcout << "Hi" << std::endl;
    
    SMCOutput* output = alg->run();

    double log_likelihood = output->log_likelihood;
    
    delete output;
    delete alg;
    
    return log_likelihood;
    
  }
  else
  {
    SMCMCMCMove* alg = new SMCMCMCMove(&rng,
                                       &seed,
                                       &the_data,
                                       algorithm_parameters,
                                       number_of_particles,
                                       smc_iterations_to_store,
                                       smc_iterations_to_store,
                                       the_mcmc,
                                       resampling_criterion,
                                       adaptive_target_criterion,
                                       sequencer_info[3],
                                       smc_termination,
                                       sequencer_variables,
                                       sequencer_schedules,
                                       likelihood_estimators,
                                       proposal_in,
                                       without_cancelled_index,
                                       full_index,
                                       proposal_is_evaluated_in,
                                       true,
                                       true,
                                       false,
                                       false,
                                       parallel_in,
                                       grain_size_in,
                                       "");
    
    //std::chrono::high_resolution_clock::time_point start_time, end_time;
    //start_time = std::chrono::high_resolution_clock::now();
    
    SMCOutput* output = alg->run();
    
    //end_time = std::chrono::high_resolution_clock::now();
    
    //std::chrono::duration<double> elapsed_time = end_time - start_time;
    
    //Rcout << elapsed_time.count() << std::endl;
    //output->set_time(elapsed_time.count());
    
    if (strcmp(results_name_in.get_cstring(),"") != 0)
      output->write(results_name_in.get_cstring());
    
    double log_likelihood = output->log_likelihood;
    
    delete output;
    delete alg;
    
    return log_likelihood;
    
  }

}

// [[Rcpp::export]]
double do_enki(const List &model,
               const List &parameters,
               const List &algorithm_parameter_list,
               size_t number_of_ensemble_members,
               const List &mcmc_termination_method,
               const List &mcmc_weights_method,
               const List &enk_sequencer_method,
               const List &adaptive_target_method,
               const List &enk_termination_method,
               const List &enk_likelihood_index_method,
               const List &enk_shifter_method,
               size_t enk_iterations_to_store,
               bool write_to_file_at_each_iteration,
               bool parallel_in,
               size_t grain_size_in,
               const String &results_name_in,
               size_t seed)
{
  RandomNumberGenerator rng;
  
  Data the_data = get_data(model);
  std::vector<Data> data_created_in_get_likelihood_estimators;
  std::vector<Data> data_created_in_get_measurement_covariance_estimators;
  
  //std::string results_name = "/Users/richard/Dropbox/code/ilike/experiments/test";
  
  // May need to alter for cases where the likelihood needs to be tuned automatically (e.g. in ABC).
  
  // Check if the prior is the proposal: affects what llhd_estimators we include.
  
  std::vector<LikelihoodEstimator*> likelihood_estimators;
  std::vector<MeasurementCovarianceEstimator*> estimators;
  IndependentProposalKernel* proposal_in;
  bool proposal_is_evaluated_in;
  
  List sequencer_info = get_smc_sequencer_info(model,
                                               parameters,
                                               enk_sequencer_method);
  std::vector<std::string> sequencer_types = sequencer_info[0];
  std::vector<std::string> sequencer_variables = sequencer_info[1];
  std::vector<std::vector<double>> sequencer_schedules = sequencer_info[2];
  
  std::vector<int> factors_affected_by_smc_sequence;
  if (sequencer_info.length()==5)
  {
    NumericVector factors_affected_by_smc_sequence_temp = sequencer_info[4];
    factors_affected_by_smc_sequence = as<std::vector<int>>(factors_affected_by_smc_sequence_temp);
  }
  
  Index* without_cancelled_index = NULL;
  Index* full_index = NULL;
  
  bool any_annealing = false;
  
  proposal_in = get_prior_as_simulate_only_proposal(model,
                                                    parameters);
  
  std::shared_ptr<Transform> transform = NULL;
  
  likelihood_estimators = get_likelihood_estimators(&rng,
                                                    &seed,
                                                    &the_data,
                                                    model,
                                                    parameters,
                                                    true,
                                                    sequencer_types,
                                                    sequencer_variables,
                                                    sequencer_schedules,
                                                    proposal_in,
                                                    without_cancelled_index,
                                                    full_index,
                                                    any_annealing,
                                                    factors_affected_by_smc_sequence,
                                                    data_created_in_get_likelihood_estimators);
  
  estimators = get_measurement_covariance_estimators(&rng,
                                                     &seed,
                                                     &the_data,
                                                     model,
                                                     parameters,
                                                     sequencer_types,
                                                     sequencer_variables,
                                                     sequencer_schedules,
                                                     enk_likelihood_index_method,
                                                     transform,
                                                     data_created_in_get_measurement_covariance_estimators);
  
  // note that this will always result in the temperature being adapted in the SMC algorithm (not smcfixed)
  /*
   likelihood_estimators = convent_to_annealed_likelihoods_if_needed(&rng,
   &seed,
   &the_data,
   sequencer_info[0],
   sequencer_info[1],
   likelihood_estimators,
   proposal_in,
   proposal_is_evaluated_in);
   */
  
  Parameters algorithm_parameters = make_algorithm_parameters(algorithm_parameter_list);
  
  MCMC* the_mcmc = make_mcmc(model,
                             parameters,
                             &the_data,
                             mcmc_termination_method,
                             mcmc_weights_method,
                             full_index);
  
  
  SMCCriterion* adaptive_target_criterion = get_adaptive_target_method(model,
                                                                       parameters,
                                                                       adaptive_target_method,
                                                                       number_of_ensemble_members);
  
  SMCTermination* enk_termination = get_smc_termination_method(model,
                                                               parameters,
                                                               enk_termination_method);
  
  EnsembleShifter* shifter = get_enk_shifter_method(model,
                                                    parameters,
                                                    enk_shifter_method);
  
  /*
   if (without_cancelled_index==NULL)
   {
   Rcout << "Without cancelled index is NULL." << std::endl;
   }
   
   if (full_index==NULL)
   {
   Rcout << "Full index is NULL." << std::endl;
   }
   */
  
  if (write_to_file_at_each_iteration)
  {
    EnsembleKalmanInversion* alg = new EnsembleKalmanInversion(&rng,
                                                               &seed,
                                                               &the_data,
                                                               number_of_ensemble_members,
                                                               enk_iterations_to_store,
                                                               shifter,
                                                               adaptive_target_criterion,
                                                               sequencer_info[3],
                                                               enk_termination,
                                                               sequencer_variables[sequencer_variables.size()-1],
                                                               sequencer_schedules[sequencer_schedules.size()-1],
                                                               proposal_in,
                                                               likelihood_estimators,
                                                               estimators,
                                                               transform,
                                                               1.0,
                                                               parallel_in,
                                                               grain_size_in,
                                                               results_name_in.get_cstring());
    
    //Rcout << "Hi" << std::endl;
    
    EnsembleKalmanOutput* output = alg->run();
    
    double log_likelihood = output->log_likelihood;
    
    delete output;
    delete alg;
    
    return log_likelihood;
    
  }
  else
  {
    EnsembleKalmanInversion* alg = new EnsembleKalmanInversion(&rng,
                                                               &seed,
                                                               &the_data,
                                                               number_of_ensemble_members,
                                                               enk_iterations_to_store,
                                                               shifter,
                                                               adaptive_target_criterion,
                                                               sequencer_info[3],
                                                               enk_termination,
                                                               sequencer_variables[sequencer_variables.size()-1],
                                                               sequencer_schedules[sequencer_schedules.size()-1],
                                                               proposal_in,
                                                               likelihood_estimators,
                                                               estimators,
                                                               transform,
                                                               1.0,
                                                               parallel_in,
                                                               grain_size_in,
                                                               "");
    
    //std::chrono::high_resolution_clock::time_point start_time, end_time;
    //start_time = std::chrono::high_resolution_clock::now();
    
    EnsembleKalmanOutput* output = alg->run();
    
    //end_time = std::chrono::high_resolution_clock::now();
    
    //std::chrono::duration<double> elapsed_time = end_time - start_time;
    
    //Rcout << elapsed_time.count() << std::endl;
    //output->set_time(elapsed_time.count());
    
    if (strcmp(results_name_in.get_cstring(),"") != 0)
      output->write(results_name_in.get_cstring());
    
    double log_likelihood = output->log_likelihood;
    
    delete output;
    delete alg;
    
    return log_likelihood;
    
  }
  
}

// [[Rcpp::export]]
double do_enkmfds(const List &model,
                        const List &parameters,
                        const List &algorithm_parameter_list,
                        size_t number_of_particles,
                  double Delta_t,
                        const List &mcmc_termination_method,
                  const List &mcmc_weights_method,
                        const List &adaptive_resampling_method,
                        const List &smc_sequencer_method,
                        const List &adaptive_target_method,
                        const List &smc_termination_method,
                        size_t smc_iterations_to_store,
                        bool write_to_file_at_each_iteration,
                        bool parallel_in,
                        size_t grain_size_in,
                        const String &results_name_in,
                        size_t seed)
{
  
  RandomNumberGenerator rng;
  
  Data the_data = get_data(model);
  std::vector<Data> data_created_in_get_likelihood_estimators;
  
  //std::string results_name = "/Users/richard/Dropbox/code/ilike/experiments/test";
  
  // May need to alter for cases where the likelihood needs to be tuned automatically (e.g. in ABC).
  
  // Check if the prior is the proposal: affects what llhd_estimators we include.
  
  std::vector<LikelihoodEstimator*> likelihood_estimators;
  IndependentProposalKernel* proposal_in;
  bool proposal_is_evaluated_in;
  
  List sequencer_info = get_smc_sequencer_info(model,
                                               parameters,
                                               smc_sequencer_method);
  std::vector<std::string> sequencer_types = sequencer_info[0];
  std::vector<std::string> sequencer_variables = sequencer_info[1];
  std::vector<std::vector<double>> sequencer_schedules = sequencer_info[2];
  
  std::vector<int> factors_affected_by_smc_sequence;
  if (sequencer_info.length()==5)
  {
    NumericVector factors_affected_by_smc_sequence_temp = sequencer_info[4];
    factors_affected_by_smc_sequence = as<std::vector<int>>(factors_affected_by_smc_sequence_temp);
  }
  
  Index* without_cancelled_index = NULL;
  Index* full_index = NULL;
  
  bool any_annealing = false;
  
  if ( model.containsElementNamed("importance_proposal") )
  {
    proposal_in = get_proposal(model,
                               parameters,
                               &the_data);
    
    likelihood_estimators = get_likelihood_estimators(&rng,
                                                      &seed,
                                                      &the_data,
                                                      model,
                                                      parameters,
                                                      true,
                                                      sequencer_types,
                                                      sequencer_variables,
                                                      sequencer_schedules,
                                                      proposal_in,
                                                      without_cancelled_index,
                                                      full_index,
                                                      any_annealing,
                                                      factors_affected_by_smc_sequence,
                                                      data_created_in_get_likelihood_estimators);
    
    proposal_is_evaluated_in = true;
    
  }
  else
  {
    Rcout << "No importance sampling proposal specified; using prior." << std::endl;
    
    proposal_in = get_prior_as_simulate_only_proposal(model,
                                                      parameters);
    
    likelihood_estimators = get_likelihood_estimators(&rng,
                                                      &seed,
                                                      &the_data,
                                                      model,
                                                      parameters,
                                                      false,
                                                      sequencer_types,
                                                      sequencer_variables,
                                                      sequencer_schedules,
                                                      proposal_in,
                                                      without_cancelled_index,
                                                      full_index,
                                                      any_annealing,
                                                      factors_affected_by_smc_sequence,
                                                      data_created_in_get_likelihood_estimators);
    
    proposal_is_evaluated_in = false;
    
  }
  
  // note that this will always result in the temperature being adapted in the SMC algorithm (not smcfixed)
  /*
   likelihood_estimators = convent_to_annealed_likelihoods_if_needed(&rng,
   &seed,
   &the_data,
   sequencer_info[0],
   sequencer_info[1],
   likelihood_estimators,
   proposal_in,
   proposal_is_evaluated_in);
   */
  
  Parameters algorithm_parameters = make_algorithm_parameters(algorithm_parameter_list);
  
  MCMC* the_mcmc = make_mcmc(model,
                             parameters,
                             &the_data,
                             mcmc_termination_method,
                             mcmc_weights_method,
                             full_index);
  
  SMCCriterion* resampling_criterion = get_resampling_method(model,
                                                             parameters,
                                                             adaptive_resampling_method,
                                                             number_of_particles);
  
  SMCCriterion* adaptive_target_criterion = get_adaptive_target_method(model,
                                                                       parameters,
                                                                       adaptive_target_method,
                                                                       number_of_particles);
  
  SMCTermination* smc_termination = get_smc_termination_method(model,
                                                               parameters,
                                                               smc_termination_method);
  
  /*
   if (without_cancelled_index==NULL)
   {
   Rcout << "Without cancelled index is NULL." << std::endl;
   }
   
   if (full_index==NULL)
   {
   Rcout << "Full index is NULL." << std::endl;
   }
   */
  
  if (write_to_file_at_each_iteration)
  {
    SMCMCMCMove* alg = new SMCMCMCMove(&rng,
                                       &seed,
                                       &the_data,
                                       algorithm_parameters,
                                       number_of_particles,
                                       smc_iterations_to_store,
                                       smc_iterations_to_store,
                                       the_mcmc,
                                       resampling_criterion,
                                       adaptive_target_criterion,
                                       sequencer_info[3],
                                       smc_termination,
                                       sequencer_variables,
                                       sequencer_schedules,
                                       likelihood_estimators,
                                       proposal_in,
                                       without_cancelled_index,
                                       full_index,
                                       proposal_is_evaluated_in,
                                       true,
                                       true,
                                       false,
                                       false,
                                       parallel_in,
                                       grain_size_in,
                                       results_name_in.get_cstring());
    
    //Rcout << "Hi" << std::endl;
    
    SMCOutput* output = alg->run();
    
    double log_likelihood = output->log_likelihood;
    
    delete output;
    delete alg;
    
    return log_likelihood;
    
  }
  else
  {
    SMCMCMCMove* alg = new SMCMCMCMove(&rng,
                                       &seed,
                                       &the_data,
                                       algorithm_parameters,
                                       number_of_particles,
                                       smc_iterations_to_store,
                                       smc_iterations_to_store,
                                       the_mcmc,
                                       resampling_criterion,
                                       adaptive_target_criterion,
                                       sequencer_info[3],
                                       smc_termination,
                                       sequencer_variables,
                                       sequencer_schedules,
                                       likelihood_estimators,
                                       proposal_in,
                                       without_cancelled_index,
                                       full_index,
                                       proposal_is_evaluated_in,
                                       true,
                                       true,
                                       false,
                                       false,
                                       parallel_in,
                                       grain_size_in,
                                       "");
    
    //std::chrono::high_resolution_clock::time_point start_time, end_time;
    //start_time = std::chrono::high_resolution_clock::now();
    
    SMCOutput* output = alg->run();
    
    //end_time = std::chrono::high_resolution_clock::now();
    
    //std::chrono::duration<double> elapsed_time = end_time - start_time;
    
    //Rcout << elapsed_time.count() << std::endl;
    //output->set_time(elapsed_time.count());
    
    if (strcmp(results_name_in.get_cstring(),"") != 0)
      output->write(results_name_in.get_cstring());
    
    double log_likelihood = output->log_likelihood;
    
    delete output;
    delete alg;
    
    return log_likelihood;
    
  }
  
}

// [[Rcpp::export]]
double do_kalman_filter(const List &model,
                        const List &parameters,
                        const List &algorithm_parameter_list,
                        const List &filtering_options_list,
                        size_t kf_iterations_to_store,
                        bool write_to_file_at_each_iteration,
                        const String &results_name_in)
{

  //Data the_data = data();
  Data the_data = get_data(model);
  
  List prior_mean_and_covariance = get_prior_mean_and_covariance(model,
                                                                 parameters);
  
  arma::colvec prior_mean_in = prior_mean_and_covariance[0];
  arma::mat prior_covariance_in = prior_mean_and_covariance[1];
  
  KalmanPredictor* predictor_in = get_kalman_predictor(model,
                                                       parameters);
  
  KalmanUpdater* updater_in = get_kalman_updater(model,
                                                 parameters);
  
  List filtering_info = get_filtering_info(parameters,
                                           filtering_options_list);
  
  
  std::string index_name_in = filtering_info[0];
  size_t first_index_in = filtering_info[1];
  size_t last_index_in = filtering_info[2];
  std::string time_name_in = filtering_info[3];
  double initial_time_in = filtering_info[4];
  std::string time_diff_name_in = filtering_info[5];
  double update_time_step_in = filtering_info[6];
  size_t predictions_per_update_in = filtering_info[7];
  std::string state_name_in = filtering_info[8];
  std::string measurement_name_in = filtering_info[9];
  
  //KalmanPredictor* predictor_in = new ExactKalmanPredictor(transition_model_A,
  //                                                         transition_model_Q);
  //KalmanUpdater* updater_in = new ExactKalmanUpdater(measurement_model_H,
  //                                                   measurement_model_R);
  
  //arma::colvec prior_mean_in = initial_mean();
  //arma::mat prior_covariance_in = initial_covariance();
  
  std::vector<std::string> measurement_names_in;
  measurement_names_in.push_back(measurement_name_in);

  Parameters algorithm_parameters = make_algorithm_parameters(algorithm_parameter_list);
  
  if (write_to_file_at_each_iteration)
  {
    KalmanFilter* alg = new KalmanFilter(&the_data,
                                         kf_iterations_to_store,
                                         state_name_in,
                                         prior_mean_in,
                                         prior_covariance_in,
                                         index_name_in,
                                         time_name_in,
                                         time_diff_name_in,
                                         measurement_names_in,
                                         first_index_in,
                                         last_index_in,
                                         predictions_per_update_in,
                                         update_time_step_in,
                                         initial_time_in,
                                         true,
                                         predictor_in,
                                         updater_in,
                                         true,
                                         results_name_in.get_cstring());
    
    KalmanFilterOutput* output = alg->run();
    
    double log_likelihood = output->log_likelihood;
    
    delete output;
    delete alg;
    
    return log_likelihood;
    
  }
  else
  {
    KalmanFilter* alg = new KalmanFilter(&the_data,
                                         kf_iterations_to_store,
                                         state_name_in,
                                         prior_mean_in,
                                         prior_covariance_in,
                                         index_name_in,
                                         time_name_in,
                                         time_diff_name_in,
                                         measurement_names_in,
                                         first_index_in,
                                         last_index_in,
                                         predictions_per_update_in,
                                         update_time_step_in,
                                         initial_time_in,
                                         true,
                                         predictor_in,
                                         updater_in,
                                         true,
                                         "");
    
    KalmanFilterOutput* output = alg->run();
    
    if (strcmp(results_name_in.get_cstring(),"") != 0)
      output->write(results_name_in.get_cstring());
    
    double log_likelihood = output->log_likelihood;
    
    delete output;
    delete alg;
    
    return log_likelihood;
    
  }
  
}
