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
#include "smc_output.h"
#include "composite_independent_proposal_kernel.h"

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

List get_single_variable_two_parameter_info(const List &model_parameters,
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
                  List info = get_single_variable_two_parameter_info(model_parameters,
                                                                     current_prior,
                                                                     type);
                  
                  std::string variable = info[0];
                  double mean = info[1];
                  double sd = info[2];
                  
                  prior_factors.push_back(new GaussianDistributionFactor(variable,
                                                                         mean,
                                                                         sd));
                  
                }
                else if (type=="lnorm")
                {
                  List info = get_single_variable_two_parameter_info(model_parameters,
                                                                     current_prior,
                                                                     type);
                  
                  std::string variable = info[0];
                  double mean = info[1];
                  double sd = info[2];
                  
                  prior_factors.push_back(new LogGaussianDistributionFactor(variable,
                                                                            mean,
                                                                            sd));
                  
                }
                else if (type=="gamma")
                {
                  List info = get_single_variable_two_parameter_info(model_parameters,
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
                List info = get_single_variable_two_parameter_info(model_parameters,
                                                                   proposal_info,
                                                                   type);
                
                std::string variable = info[0];
                double mean = info[1];
                double sd = info[2];
                
                proposal = new GaussianIndependentProposalKernel(variable,
                                                                 mean,
                                                                 sd);
                
              }
              else if (type=="lnorm")
              {
                List info = get_single_variable_two_parameter_info(model_parameters,
                                                                   proposal_info,
                                                                   type);
                
                std::string variable = info[0];
                double mean = info[1];
                double sd = info[2];
                
                proposal = new LogGaussianIndependentProposalKernel(variable,
                                                                    mean,
                                                                    sd);
                
              }
              else if (type=="gamma")
              {
                List info = get_single_variable_two_parameter_info(model_parameters,
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
            proposal = new CustomGuidedDistributionProposalKernel(load_simulate_guided_distribution(simulate_importance_proposal_SEXP),
                                                            load_evaluate_log_guided_distribution(evaluate_log_importance_proposal_SEXP),data);
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
                List info = get_single_variable_two_parameter_info(model_parameters,
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
              else if (type=="lnorm")
              {
                List info = get_single_variable_two_parameter_info(model_parameters,
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
              else if (type=="gamma")
              {
                List info = get_single_variable_two_parameter_info(model_parameters,
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
  
  //std::string results_name = "/Users/richard/Dropbox/code/ilike/experiments/test";
  
  // May need to alter for cases where the likelihood needs to be tuned automatically (e.g. in ABC).
  
  // Check if the prior is the proposal: affects what llhd_estimators we include.
  
  std::vector<LikelihoodEstimator*> likelihood_estimators;
  IndependentProposalKernel* proposal_in;
  bool proposal_is_evaluated_in;
  
  if ( model.containsElementNamed("importance_proposal") )
  {
    
    likelihood_estimators = get_likelihood_eatimators(&rng,
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
    likelihood_estimators = get_likelihood_eatimators(&rng,
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
    
    likelihood_estimators = get_likelihood_eatimators(&rng,
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
