#include "reinforce_gradient_estimator.h"
#include "utils.h"
#include "gradient_estimator_output.h"
#include "reinforce_gradient_estimator_output.h"
#include "transform.h"

namespace ilike
{
ReinforceGradientEstimator::ReinforceGradientEstimator()
:GradientEstimator()
{
}

ReinforceGradientEstimator::~ReinforceGradientEstimator()
{
  
}

ReinforceGradientEstimator::ReinforceGradientEstimator(DataSubsampler* subsampler_in)
:GradientEstimator()
{
  this->subsampler = subsampler_in;
}

//Copy constructor for the ReinforceGradientEstimator class.
ReinforceGradientEstimator::ReinforceGradientEstimator(const ReinforceGradientEstimator &another)
:GradientEstimator(another)
{
  this->make_copy(another);
}

void ReinforceGradientEstimator::operator=(const ReinforceGradientEstimator &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  GradientEstimator::operator=(another);
  this->make_copy(another);
}

GradientEstimator* ReinforceGradientEstimator::duplicate() const
{
  return( new ReinforceGradientEstimator(*this));
}

void ReinforceGradientEstimator::make_copy(const ReinforceGradientEstimator &another)
{
  this->gaussian_proposal = another.gaussian_proposal;
  this->num_points = another.num_points;
  this->subsampler = another.subsampler;
  this->size_of_subsample = another.size_of_subsample;
}

GradientEstimatorOutput* ReinforceGradientEstimator::initialise()
{
  return new ReinforceGradientEstimatorOutput(this);
}

boost::unordered_map< std::string, std::vector<arma::mat>> ReinforceGradientEstimator::simulate_auxiliary_variables()
{
  boost::unordered_map< std::string, std::vector<arma::mat>> output;
  
  std::vector<std::string> variables = this->gaussian_proposal.get_variables();
  for (auto j=variables.begin(); j!=variables.end(); ++j)
  {
    std::vector<arma::mat> offsets;
    offsets.reserve(this->num_points);
    output[*j] = offsets;
  }
  
  for (size_t i=0; i<this->num_points; ++i)
  {
    Parameters simulated_offsets = this->gaussian_proposal.independent_simulate(*this->rng);
    
    for (auto j=variables.begin(); j!=variables.end(); ++j)
    {
      output[*j].push_back(simulated_offsets[*j]);
    }
  }
  
  return output;
}

arma::mat ReinforceGradientEstimator::get_gradient_of_log(const std::string &variable,
                                                          const std::vector<arma::mat> &auxiliary_variables,
                                                          const Index* index,
                                                          const Particle &particle)
{
  // could change to use mean of points sampled from proposal
  Particle current_particle = particle;
  arma::mat current_point = current_particle.parameters[variable];
  
  this->subsampler->subsample(this->size_of_subsample);
  
  double current_subsampled_target = current_particle.subsample_evaluate_likelihoods(index);
  
  double jacobian_determinant_of_particle;
  
  if (this->proposal->transform!=NULL)
  {
    jacobian_determinant_of_particle = this->proposal->transform->log_abs_inverse_jacobian_determinant(current_particle.get_transformed_parameters(this->proposal));
  }
  else
  {
    jacobian_determinant_of_particle = 1.0;
  }
  
  arma::mat output(current_point.n_rows,current_point.n_cols);
  
  for (auto i=auxiliary_variables.begin();
       i!=auxiliary_variables.end();
       ++i)
  {
    arma::mat simulated_point = current_point + *i;
    Particle simulated_particle = particle;
    simulated_particle.parameters[variable] = simulated_point;
    simulated_particle.simulate_factor_variables();
    simulated_particle.simulate_ensemble_factor_variables();
    
    double jacobian_determinant_of_simulated;
    
    if (this->proposal->transform!=NULL)
    {
      // project simulated particle back into original space (change params and move_params, I think)
      
      // multiply result by Jacobian det of inverse
      jacobian_determinant_of_simulated = this->proposal->transform->log_abs_inverse_jacobian_determinant(simulated_particle.get_transformed_parameters(this->proposal));
    }
    else
    {
      jacobian_determinant_of_simulated = 1.0;
    }
    
    output = output + (jacobian_determinant_of_simulated + simulated_particle.subsample_evaluate_likelihoods(index) - jacobian_determinant_of_particle - current_particle.subsample_target_evaluated)*(this->gaussian_proposal.get_inverse_covariance(variable)*(*i));
  }
  return (1.0/double(this->num_points))*output;
}

/*
 arma::mat ReinforceGradientEstimator::get_gradient_of_log(const std::string &variable,
 const Index* index,
 Particle &particle,
 const Parameters &conditioned_on_parameters)
 {
 // could change to use mean of points sampled from proposal
 this->subsampler->subsample(this->size_of_subsample);
 
 //arma::colvec unnormalised_weights(this->num_points);
 //std::vector<arma::mat> contributions;
 //contributions.reserve(this->num_points);
 
 Particle simulated_particle = this->gaussian_proposal.subsample_move(*this->rng,
 variable,
 particle,
 conditioned_on_parameters);
 
 arma::mat output;
 double jacobian_determinant_of_simulated;
 double jacobian_determinant_of_particle;
 if (this->proposal->transform!=NULL)
 {
 // project simulated particle back into original space (change params and move_params, I think)
 simulated_particle.move_transformed_parameters = simulated_particle.parameters;
 simulated_particle.parameters = this->proposal->transform->inverse_transform(simulated_particle.move_transformed_parameters);
 // do evaluation
 // multiply result by Jacobian det of inverse
 jacobian_determinant_of_simulated = this->proposal->transform->log_abs_inverse_jacobian_determinant(simulated_particle.move_transformed_parameters);
 // particle.subsample_target_evaluated also needs to be multiplied by Jacobian det at its value (adding in log space)
 jacobian_determinant_of_particle = this->proposal->transform->log_abs_inverse_jacobian_determinant(particle.move_transformed_parameters);
 }
 else
 {
 jacobian_determinant_of_simulated = 1.0;
 jacobian_determinant_of_particle = 1.0;
 }
 
 output = (jacobian_determinant_of_simulated + simulated_particle.subsample_evaluate_likelihoods(index,conditioned_on_parameters) - jacobian_determinant_of_particle - particle.subsample_target_evaluated)*(this->gaussian_proposal.get_inverse_covariance(variable)*(particle.parameters[variable]-simulated_particle.parameters[variable]));
 
 for (size_t i=1; i<this->num_points; ++i)
 {
 simulated_particle = this->gaussian_proposal.subsample_move(*this->rng,
 variable,
 particle,
 conditioned_on_parameters);
 
 if (this->proposal->transform!=NULL)
 {
 // project simulated particle back into original space (change params and move_params, I think)
 simulated_particle.move_transformed_parameters = simulated_particle.parameters;
 simulated_particle.parameters = this->proposal->transform->inverse_transform(simulated_particle.move_transformed_parameters);
 // do evaluation
 // multiply result by Jacobian det of inverse
 jacobian_determinant_of_simulated = this->proposal->transform->log_abs_inverse_jacobian_determinant(simulated_particle.move_transformed_parameters);
 }
 else
 {
 jacobian_determinant_of_simulated = 1.0;
 }
 
 output = (jacobian_determinant_of_simulated + simulated_particle.subsample_evaluate_likelihoods(index,conditioned_on_parameters) - jacobian_determinant_of_particle - particle.subsample_target_evaluated)*(this->gaussian_proposal.get_inverse_covariance(variable)*(particle.parameters[variable]-simulated_particle.parameters[variable]));
 }
 return (1.0/double(this->num_points))*output;
 }
 */

arma::mat ReinforceGradientEstimator::subsample_get_gradient_of_log(const std::string &variable,
                                                                    const std::vector<arma::mat> &auxiliary_variables,
                                                                    const Index* index,
                                                                    const Particle &particle)
{
  // not used, since would need to have subsample of subsample, and need a different approach to doing the scaling
  Rcpp::stop("ReinforceGradientEstimator::subsample_get_gradient_of_log - not implemented.");
}

/*
 arma::mat ReinforceGradientEstimator::subsample_get_gradient_of_log(const std::string &variable,
 const Index* index,
 Particle &particle,
 const Parameters &conditioned_on_parameters)
 {
 // not used, since would need to have subsample of subsample, and need a different approach to doing the scaling
 Rcpp::stop("ReinforceGradientEstimator::subsample_get_gradient_of_log - not implemented.");
 }
 */
}
