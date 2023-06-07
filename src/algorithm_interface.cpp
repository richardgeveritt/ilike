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
    return load_data(data_SEXP);
  }
  else
  {
    Rcpp::stop("do_importance_sampler: data not found in model specification.");
  }
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
      Rcpp::stop("Parameter index does not correspond to a real number in the parameters file.");
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
      Rcpp::stop("Parameter index does not correspond to a real number in the parameters file.");
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
      Rcpp::stop("Parameter index does not correspond to a real number in the parameters file.");
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

std::vector<LikelihoodEstimator*> get_likelihood_estimators(RandomNumberGenerator* rng_in,
                                                            size_t* seed_in,
                                                            Data* data_in,
                                                            const List &model,
                                                            const List &model_parameters,
                                                            bool include_priors)
{
  
  std::vector<LikelihoodEstimator*> likelihood_estimators;
  std::vector<DistributionFactor*> numerator_distribution_factors;
  std::vector<LikelihoodFactor*> numerator_likelihood_factors;
  
  std::vector<DistributionFactor*> prior_factors;
  std::vector<LikelihoodFactor*> exact_likelihood_factors;

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
          if (include_priors)
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
                  
                  prior_factors.push_back(new GaussianDistributionFactor(variable,
                                                                         mean,
                                                                         sd));
                  
                }
                else if (type=="mvnorm")
                {
                  List info = get_single_variable_vector_and_matrix_parameter_info(model_parameters,
                                                                                   current_prior,
                                                                                   type);
                  
                  std::string variable = info[0];
                  arma::colvec mean = info[1];
                  arma::mat cov = info[2];
                  
                  prior_factors.push_back(new GaussianDistributionFactor(variable,
                                                                         mean,
                                                                         cov));
                  
                }
                else if (type=="lnorm")
                {
                  List info = get_single_variable_two_double_parameter_info(model_parameters,
                                                                            current_prior,
                                                                            type);
                  
                  std::string variable = info[0];
                  double mean = info[1];
                  double sd = info[2];
                  
                  prior_factors.push_back(new LogGaussianDistributionFactor(variable,
                                                                            mean,
                                                                            sd));
                  
                }
                else if (type=="mvlnorm")
                {
                  List info = get_single_variable_vector_and_matrix_parameter_info(model_parameters,
                                                                                   current_prior,
                                                                                   type);
                  
                  std::string variable = info[0];
                  arma::colvec mean = info[1];
                  arma::mat cov = info[2];
                  
                  prior_factors.push_back(new LogGaussianDistributionFactor(variable,
                                                                            mean,
                                                                            cov));
                  
                }
                else if (type=="gamma")
                {
                  List info = get_single_variable_two_double_parameter_info(model_parameters,
                                                                     current_prior,
                                                                     type);
                  
                  std::string variable = info[0];
                  double shape = info[1];
                  double rate = info[2];
                  
                  prior_factors.push_back(new GammaDistributionFactor(variable,
                                                                      shape,
                                                                      rate));
                  
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
        }
        else if (current_factor.containsElementNamed("evaluate_log_prior"))
        {
          if (include_priors)
          {
            
            SEXP evaluate_log_prior_SEXP = current_factor["evaluate_log_prior"];
            prior_factors.push_back(new CustomDistributionFactor(load_evaluate_log_distribution(evaluate_log_prior_SEXP)));
            
          }
        }
        else if ( current_factor.containsElementNamed("evaluate_log_likelihood") )
        {
          
          SEXP evaluate_log_likelihood_SEXP = current_factor["evaluate_log_likelihood"];
          exact_likelihood_factors.push_back(new CustomLikelihoodFactor(load_evaluate_log_likelihood(evaluate_log_likelihood_SEXP),
                                                                        data_in));
          
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

  
  likelihood_estimators.push_back(new ExactLikelihoodEstimator(rng_in,
                                                               seed_in,
                                                               data_in,
                                                               prior_factors,
                                                               exact_likelihood_factors,
                                                               true));
  
  return likelihood_estimators;
  
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
                                                cov);
          
        }
        else if (type=="hmc")
        {
          List info = get_single_variable_matrix_parameter_info(model_parameters,
                                                                proposal_info,
                                                                type);
          
          std::string variable = info[0];
          arma::mat cov = info[1];
          
          proposal = new HMCProposalKernel(variable,
                                           cov);
          
        }
        else if (type=="barker_dynamics")
        {
          List info = get_single_variable_matrix_parameter_info(model_parameters,
                                                                proposal_info,
                                                                type);
          
          std::string variable = info[0];
          arma::mat cov = info[1];
          
          proposal = new BarkerDynamicsProposalKernel(variable,
                                                      cov);
          
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

MCMC* make_mcmc(const List &model,
                const List &model_parameters,
                const Data* data,
                size_t number_of_mcmc_iterations)
{
  MCMC* mcmc = NULL;
  
  if (model.containsElementNamed("order_of_mcmc"))
  {
    NumericVector order_of_mcmc = model["order_of_mcmc"];
    std::vector<MCMC*> moves;
    moves.reserve(order_of_mcmc.size());
    
    size_t simulate_mh_proposal_index = 0;
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
          
          mcmc = new MetropolisHastingsMCMC(number_of_mcmc_iterations,
                                            proposal);
          
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
          
          moves.push_back(new MetropolisHastingsMCMC(number_of_mcmc_iterations,
                                                     proposal));
          
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
          
          moves.push_back(new MetropolisMCMC(number_of_mcmc_iterations,
                                             proposal));
          
          simulate_m_proposal_index = simulate_m_proposal_index + 1;
        }
        else
        {
          stop("No m_proposal found.");
        }
        
      }
      else
      {
        Rcpp::stop("get_mcmc: invalid type in order_of_mcmc.");
      }
      
      moves.push_back(mcmc);
      
    }
    
    if (model.containsElementNamed("mcmc_weights"))
    {
      // do stochastic if more than one
      
      NumericVector mcmc_weights;
      List mcmc_weights_list = model["mcmc_weights"];
      List mcmc_weights_list0 = mcmc_weights_list[0];
      SEXP mcmc_weights_SEXP = mcmc_weights_list0["mcmc_weights"];
      mcmc_weights = load_mcmc_weights(mcmc_weights_SEXP);
      
      if (moves.size()==1)
      {
        Rcout << "one move" << std::endl;
        return mcmc;
      }
      else
      {
        Rcout << "multiple moves" << std::endl;
        return new StochasticScanMCMC(moves,
                                      mcmc_weights);
      }
      
    }
    else
    {
      if (moves.size()==1)
      {
        return mcmc;
      }
      else
      {
        return new DeterministicScanMCMC(moves);
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
    List list_params = initial_values[i];
    if (list_params.size()>0)
    {
      CharacterVector names = list_params.names();
      
      for (size_t i=0; i<names.size(); ++i)
      {
        arma::mat param = list_params[Rcpp::as<std::string>(names[i])];
        current_parameters[Rcpp::as<std::string>(names[i])] = param;
      }
    }
    else
    {
      Rcpp::stop("List of initial values has empty entries.");
    }
    
    output.push_back(current_parameters);
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
  return output;
}

// [[Rcpp::export]]
void do_importance_sampler(const List &model,
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
  
  //std::string results_name = "/Users/richard/Dropbox/code/ilike/experiments/test";
  
  // May need to alter for cases where the likelihood needs to be tuned automatically (e.g. in ABC).
  
  // Check if the prior is the proposal: affects what llhd_estimators we include.
  
  std::vector<LikelihoodEstimator*> likelihood_estimators;
  IndependentProposalKernel* proposal_in;
  bool proposal_is_evaluated_in;
  
  if ( model.containsElementNamed("importance_proposal") )
  {
    
    likelihood_estimators = get_likelihood_estimators(&rng,
                                                      &seed,
                                                      &the_data,
                                                      model,
                                                      parameters,
                                                      true);
    
    proposal_in = get_proposal(model,
                               parameters,
                               &the_data);
    
    proposal_is_evaluated_in = true;
    
  }
  else
  {
    Rcout << "No importance sampling proposal specified; using prior." << std::endl;
    
    likelihood_estimators = get_likelihood_estimators(&rng,
                                                      &seed,
                                                      &the_data,
                                                      model,
                                                      parameters,
                                                      false);
    
    proposal_in = get_prior_as_simulate_only_proposal(model,
                                                      parameters);
    
    proposal_is_evaluated_in = false;
  }
  
  
  Parameters algorithm_parameters = make_algorithm_parameters(algorithm_parameter_list);
  
  ImportanceSampler alg(&rng,
                        &seed,
                        &the_data,
                        algorithm_parameters,
                        number_of_importance_points,
                        likelihood_estimators,
                        proposal_in,
                        proposal_is_evaluated_in,
                        true,
                        true,
                        parallel_in,
                        grain_size_in,
                        "");
  
  std::chrono::high_resolution_clock::time_point start_time, end_time;
  start_time = std::chrono::high_resolution_clock::now();
  
  SMCOutput* output = alg.run();
  
  end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_time = end_time - start_time;
  output->set_time(elapsed_time.count());
  
  if (strcmp(results_name_in.get_cstring(),"") != 0)
    output->write(results_name_in.get_cstring());
  delete output;
  
}
  
// [[Rcpp::export]]
void do_mcmc(const List &model,
             const List &parameters,
             const List &algorithm_parameter_list,
             const List &initial_values,
             size_t number_of_mcmc_iterations,
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
  
  // May need to alter for cases where the likelihood needs to be tuned automatically (e.g. in ABC).
  
  // Check if the prior is the proposal: affects what llhd_estimators we include.
  
  
  std::vector<LikelihoodEstimator*> likelihood_estimators;
  likelihood_estimators = get_likelihood_estimators(&rng,
                                                    &seed,
                                                    &the_data,
                                                    model,
                                                    parameters,
                                                    true);
  
  Parameters algorithm_parameters = make_algorithm_parameters(algorithm_parameter_list);
  
  
  
  MCMC* the_mcmc = make_mcmc(model,
                             parameters,
                             &the_data,
                             number_of_mcmc_iterations);
  
  
  SMCMCMCMove* alg;
  
  if (initial_values.size()==0)
  {
    
    IndependentProposalKernel* proposal_in = get_prior_as_proposal(model,
                                                                   parameters);
    
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
                          parallel_in,
                          grain_size_in,
                          "");
    
  }
  
  
  std::chrono::high_resolution_clock::time_point start_time, end_time;
  start_time = std::chrono::high_resolution_clock::now();
  
  SMCOutput* output = alg->run();
  
  end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_time = end_time - start_time;
  output->set_time(elapsed_time.count());
  
  if (strcmp(results_name_in.get_cstring(),"") != 0)
    output->write(results_name_in.get_cstring());
  delete output;
  
  delete alg;
}


// [[Rcpp::export]]
void do_smc_mcmc(const List &model,
                 const List &parameters,
                 size_t number_of_importance_points,
                 bool parallel_in,
                 size_t grain_size_in,
                 const String &results_name_in,
                 size_t seed)
{
  
  RandomNumberGenerator rng;
  
  Data the_data = get_data(model);
  
  //std::string results_name = "/Users/richard/Dropbox/code/ilike/experiments/test";
  
  // May need to alter for cases where the likelihood needs to be tuned automatically (e.g. in ABC).
  
  // Check if the prior is the proposal: affects what llhd_estimators we include.
  
  
  std::vector<LikelihoodEstimator*> likelihood_estimators;
  IndependentProposalKernel* proposal_in;
  bool proposal_is_evaluated_in;
  
  if ( model.containsElementNamed("importance_proposals") )
  {
    likelihood_estimators = get_likelihood_estimators(&rng,
                                                      &seed,
                                                      &the_data,
                                                      model,
                                                      parameters,
                                                      true);
    
    proposal_in = get_proposal(model,
                               parameters,
                               &the_data);
    
    proposal_is_evaluated_in = true;
  }
  else
  {
    Rcout << "No importance sampling proposal specified; using prior." << std::endl;
    
    likelihood_estimators = get_likelihood_estimators(&rng,
                                                      &seed,
                                                      &the_data,
                                                      model,
                                                      parameters,
                                                      false);
    /*
    proposal_in = get_prior_as_simulate_only_proposal(model,
                                                      parameters);
    
    proposal_is_evaluated_in = false;
    */
  }
  
  /*
  ImportanceSampler alg(&rng,
                        &seed,
                        &the_data,
                        number_of_importance_points,
                        likelihood_estimators,
                        proposal_in,
                        proposal_is_evaluated_in,
                        true,
                        true,
                        parallel_in,
                        grain_size_in,
                        "");
  
  std::chrono::high_resolution_clock::time_point start_time, end_time;
  start_time = std::chrono::high_resolution_clock::now();
  
  SMCOutput* output = alg.run();
  
  end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_time = end_time - start_time;
  output->set_time(elapsed_time.count());
  
  if (strcmp(results_name_in.get_cstring(),"") != 0)
    output->write(results_name_in.get_cstring());
  delete output;
  
  Rcout << "Sampling finished." << std::endl;
  if (strcmp(results_name_in.get_cstring(),"") != 0)
    Rcout << "Results written to " << results_name_in.get_cstring() << "_smc" << std::endl;
  */
}
