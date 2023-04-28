#include <RcppArmadillo.h>
using namespace Rcpp;

#include <string>
#include <sstream>
#include <vector>
#include <chrono>

#include "distributions.h"
#include "parameters.h"
#include "exact_likelihood_estimator.h"
#include "importance_sampler.h"
#include "gaussian_distribution_factor.h"
#include "custom_distribution_factor.h"
#include "custom_likelihood_factor.h"
#include "gaussian_independent_proposal_kernel.h"
#include "custom_distribution_proposal_kernel.h"
#include "smc_output.h"

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
  if (index_from_one>list.size())
  {
    SEXP var = list[index_from_one-1];
    return Rf_isReal(var) && Rf_length(var) == 1;
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
    SEXP data_SEXP = data_list[0];
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

std::vector<LikelihoodEstimator*> get_likelihood_eatimators(RandomNumberGenerator* rng_in,
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
  
  if (include_priors)
  {
    if ( model.containsElementNamed("prior") )
    {
      List priors = model["prior"];
      
      for (size_t i=0; i<priors.size(); ++i)
      {
        if (Rf_isNewList(priors[i]))
        {
          List current_prior = priors[i];
          if ( current_prior.containsElementNamed("type") && current_prior.containsElementNamed("variables") )
          {
            std::string type = current_prior["type"];
            if (type=="norm")
            {
              std::string augmented_variable_names = current_prior["variables"];
              
              // split string in ;
              std::vector<std::string> variable_names = split(augmented_variable_names,';');
              
              // throw error if more than one variable
              if (variable_names.size()!=1)
                Rcpp::stop("Only one variable allowed for univariate Gaussian distribution.");
              
              if (!current_prior.containsElementNamed("parameters"))
                stop("Parameters missing for univariate Gaussian prior.");
              
              List parameters = current_prior["parameters"];
              
              if (parameters.size()!=2)
                Rcpp::stop("Mean and standard deviation parameters required for univariate Gaussian prior.");
              
              double mean = extract_double_parameter(parameters,
                                                     model_parameters,
                                                     0);
              
              double sd = extract_double_parameter(parameters,
                                                   model_parameters,
                                                   1);
              
              prior_factors.push_back(new GaussianDistributionFactor(variable_names[0],
                                                                     mean,
                                                                     sd));
              
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
          stop("Error in prior section of model file.");
        }
      }
    }
    
    if ( model.containsElementNamed("log_evaluate_log_prior") )
    {
      List priors = model["evaluate_log_prior"];
      
      for (size_t i=0; i<priors.size(); ++i)
      {
        SEXP prior_log_evaluate_SEXP = priors[i];
        prior_factors.push_back(new CustomDistributionFactor(load_evaluate_log_distribution(prior_log_evaluate_SEXP)));
      }
    }
  }
  
  if ( model.containsElementNamed("evaluate_log_likelihood") )
  {
    List exact_likelihoods = model["evaluate_log_likelihood"];
    
    for (size_t i=0; i<exact_likelihoods.size(); ++i)
    {
      SEXP likelihood_log_evaluate_SEXP = exact_likelihoods[i];
      exact_likelihood_factors.push_back(new CustomLikelihoodFactor(load_evaluate_log_likelihood(likelihood_log_evaluate_SEXP),
                                                                    data_in));
    }
    
  }

  likelihood_estimators.push_back(new ExactLikelihoodEstimator(rng_in,
                                                               seed_in,
                                                               data_in,
                                                               prior_factors,
                                                               exact_likelihood_factors,
                                                               true));

  return likelihood_estimators;
}

IndependentProposalKernel* get_proposal(const List &model,
                                        const List &model_parameters)
{
  IndependentProposalKernel* proposal = NULL;
  
  if ( model.containsElementNamed("proposal") && model.containsElementNamed("evaluate_log_proposal") )
  {
    stop("Model file contains a 'proposal' section and a 'evaluate_log_proposal' section. Include only one of these.");
  }
  
  if ( model.containsElementNamed("proposal") && model.containsElementNamed("proposal_simulate") )
  {
    stop("Model file contains a 'proposal' section and a 'simulate_proposal' section. Include only one of these.");
  }
  
  if ( model.containsElementNamed("proposal") )
  {
    List proposal_info = model["proposal"];
    
    if (Rf_isNewList(proposal_info))
    {
      if ( proposal_info.containsElementNamed("type") && proposal_info.containsElementNamed("variables") )
      {
        std::string type = proposal_info["type"];
        if (type=="norm")
        {
          std::string augmented_variable_names = proposal_info["variables"];
          
          // split string in ;
          std::vector<std::string> variable_names = split(augmented_variable_names,';');
          
          // throw error if more than one variable
          if (variable_names.size()!=1)
            Rcpp::stop("Only one variable allowed for univariate Gaussian distribution.");
          
          if (!proposal_info.containsElementNamed("parameters"))
            stop("Parameters missing for univariate Gaussian proposal.");
          
          List parameters = proposal_info["parameters"];
          
          if (parameters.size()!=2)
            Rcpp::stop("Mean and standard deviation parameters required for univariate Gaussian proposal.");
          
          double mean = extract_double_parameter(parameters,
                                                 model_parameters,
                                                 0);
          
          double sd = extract_double_parameter(parameters,
                                               model_parameters,
                                               1);
          
          proposal = new GaussianIndependentProposalKernel(variable_names[0],
                                                           mean,
                                                           sd);
          
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
  else if ( model.containsElementNamed("evaluate_log_proposal") && model.containsElementNamed("simulate_proposal") )
  {
    SEXP evaluate_log_proposal_SEXP = model["evaluate_log_proposal"];
    SEXP simulate_proposal_SEXP = model["simulate_proposal"];
    proposal = new CustomDistributionProposalKernel(load_simulate_distribution(simulate_proposal_SEXP),
                                                    load_evaluate_log_distribution(evaluate_log_proposal_SEXP));
  }
  else
  {
    stop("No proposal specified.");
  }
  
  return proposal;
}

IndependentProposalKernel* get_prior_as_simulate_only_proposal(const List &model,
                                                               const List &model_parameters)
{
  IndependentProposalKernel* proposal = NULL;
  
  if ( model.containsElementNamed("prior") && model.containsElementNamed("simulate_prior") )
  {
    stop("Model file contains a 'prior' section and a 'simulate_prior' section. No proposal specified, so using prior as proposal, but since both sections are specified there is ambiguity over which to use.");
  }
  
  if ( model.containsElementNamed("prior") )
  {
    List priors = model["prior"];
    
    if (priors.size()!=1)
    {
      stop("Need only a single element in the 'prior' section in the model file. No proposal specified, so using prior as proposal, but need exactly one to be specified or there is ambiguity over which to use.");
    }
    
    if (Rf_isNewList(priors[0]))
    {
      List current_prior = priors[0];
      if ( current_prior.containsElementNamed("type") && current_prior.containsElementNamed("variables") )
      {
        std::string type = current_prior["type"];
        if (type=="norm")
        {
          std::string augmented_variable_names = current_prior["variables"];
          
          // split string in ;
          std::vector<std::string> variable_names = split(augmented_variable_names,';');
          
          // throw error if more than one variable
          if (variable_names.size()!=1)
            Rcpp::stop("Only one variable allowed for univariate Gaussian distribution.");
          
          if (!current_prior.containsElementNamed("parameters"))
            stop("Parameters missing for univariate Gaussian prior.");
          
          List parameters = current_prior["parameters"];
          
          if (parameters.size()!=2)
            Rcpp::stop("Mean and standard deviation parameters required for univariate Gaussian prior.");
          
          double mean = extract_double_parameter(parameters,
                                                 model_parameters,
                                                 0);
          
          double sd = extract_double_parameter(parameters,
                                               model_parameters,
                                               1);
          
          proposal = new GaussianIndependentProposalKernel(variable_names[0],
                                                           mean,
                                                           sd);
          
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
      stop("Error in prior section of model file.");
    }
  }
  
  if ( model.containsElementNamed("simulate_prior") )
  {
    List priors = model["simulate_prior"];
    
    if (priors.size()!=1)
    {
      stop("Need only a single element in the 'simulate_prior' section in the model file. No proposal specified, so using prior as proposal, but need exactly one to be specified or there is ambiguity over which to use.");
    }
    
    SEXP simulate_prior_SEXP = priors[0];
    proposal = new CustomDistributionProposalKernel(load_simulate_distribution(simulate_prior_SEXP));
  }
  
  if (proposal==NULL)
  {
    stop("No suitable prior to use as proposal.");
  }
  
  return proposal;
}

// [[Rcpp::export]]
size_t ilike_rdtsc()
{
  return rdtsc();
}

// [[Rcpp::export]]
void do_importance_sampler(const List &model,
                           const List &parameters,
                           size_t number_of_importance_points,
                           bool parallel_in,
                           size_t grain_size_in,
                           const String &results_name_in,
                           size_t seed)
{
  
  RandomNumberGenerator rng;
  
  Data the_data = get_data(model);
  
  bool parallel = FALSE;
  bool smcfixed_flag = TRUE;
  size_t grain_size = 1;
  
  //std::string results_name = "/Users/richard/Dropbox/code/ilike/experiments/test";
  
  // May need to alter for cases where the likelihood needs to be tuned automatically (e.g. in ABC).
  
  // Check if the prior is the proposal: affects what llhd_estimators we include.
  
  std::vector<LikelihoodEstimator*> likelihood_estimators;
  IndependentProposalKernel* proposal_in;
  bool proposal_is_evaluated_in;
  
  if ( model.containsElementNamed("proposal") || model.containsElementNamed("simulate_proposal") )
  {
    likelihood_estimators = get_likelihood_eatimators(&rng,
                                                      &seed,
                                                      &the_data,
                                                      model,
                                                      parameters,
                                                      true);
    
    proposal_in = get_proposal(model,
                               parameters);
    
    proposal_is_evaluated_in = true;
  }
  else
  {
    Rcout << "No importance sampling proposal specified; using prior." << std::endl;
    
    likelihood_estimators = get_likelihood_eatimators(&rng,
                                                      &seed,
                                                      &the_data,
                                                      model,
                                                      parameters,
                                                      false);
    
    proposal_in = get_prior_as_simulate_only_proposal(model,
                                                      parameters);
    
    proposal_is_evaluated_in = false;
  }
  
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
  
  clock_t start, end;
  int max = 0;
  start = clock();
  SMCOutput* output = alg.run();
  end = clock();
  double time = double(end - start)/CLOCKS_PER_SEC;
  output->set_time(time);
  if (strcmp(results_name_in.get_cstring(),"") != 0)
    output->write(results_name_in.get_cstring());
  delete output;
  
}
