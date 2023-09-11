#include "hmm_ensemble_factor_variables.h"
#include "hmm_ensemble_factors.h"
#include "measurement_covariance_estimator_output.h"

HMMEnsembleFactorVariables::HMMEnsembleFactorVariables()
  :EnsembleFactorVariables()
{
  this->measurement_covariance_estimator_outputs.resize(0);
  this->hmm_ensemble_factors = NULL;
}

HMMEnsembleFactorVariables::~HMMEnsembleFactorVariables()
{
  for (std::vector<MeasurementCovarianceEstimatorOutput*>::iterator i=this->measurement_covariance_estimator_outputs.begin();
       i!=this->measurement_covariance_estimator_outputs.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
}

HMMEnsembleFactorVariables::HMMEnsembleFactorVariables(const HMMEnsembleFactors* hmm_ensemble_factors_in,
                                                       const std::vector<MeasurementCovarianceEstimatorOutput*> &measurement_covariance_estimator_outputs_in)
: EnsembleFactorVariables()
{
  this->measurement_covariance_estimator_outputs = measurement_covariance_estimator_outputs_in;
  this->hmm_ensemble_factors = hmm_ensemble_factors_in;
}

//Copy constructor for the HMMEnsembleFactorVariables class.
HMMEnsembleFactorVariables::HMMEnsembleFactorVariables(const HMMEnsembleFactorVariables &another)
  :EnsembleFactorVariables(another)
{
  this->make_copy(another);
}

void HMMEnsembleFactorVariables::operator=(const HMMEnsembleFactorVariables &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  for (std::vector<MeasurementCovarianceEstimatorOutput*>::iterator i=this->measurement_covariance_estimator_outputs.begin();
       i!=this->measurement_covariance_estimator_outputs.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  this->measurement_covariance_estimator_outputs.clear();
  
  EnsembleFactorVariables::operator=(another);
  this->make_copy(another);
}

EnsembleFactorVariables* HMMEnsembleFactorVariables::duplicate() const
{
  return( new HMMEnsembleFactorVariables(*this));
}

void HMMEnsembleFactorVariables::make_copy(const HMMEnsembleFactorVariables &another)
{
  this->measurement_covariance_estimator_outputs.resize(0);
  this->measurement_covariance_estimator_outputs.reserve(another.measurement_covariance_estimator_outputs.size());
  for (std::vector<MeasurementCovarianceEstimatorOutput*>::const_iterator i=another.measurement_covariance_estimator_outputs.begin();
       i!=another.measurement_covariance_estimator_outputs.end();
       ++i)
  {
    if (*i!=NULL)
      this->measurement_covariance_estimator_outputs.push_back((*i)->duplicate());
    else
      this->measurement_covariance_estimator_outputs.push_back(NULL);
  }
  this->hmm_ensemble_factors = another.hmm_ensemble_factors;
}

std::vector<arma::rowvec> HMMEnsembleFactorVariables::get_measurement_states_for_covariance() const
{
  std::vector<arma::rowvec> measurements;
  measurements.reserve(this->measurement_covariance_estimator_outputs.size());
  for (auto i=this->measurement_covariance_estimator_outputs.begin();
       i!=this->measurement_covariance_estimator_outputs.end();
       ++i)
  {
    measurements.push_back((*i)->get_measurement_state_for_covariance());
  }
  return measurements;
}

std::vector<arma::colvec> HMMEnsembleFactorVariables::get_shifts(double inverse_incremental_temperature) const
{
  std::vector<arma::colvec> shifts;
  shifts.reserve(this->measurement_covariance_estimator_outputs.size());
  for (auto i=this->measurement_covariance_estimator_outputs.begin();
       i!=this->measurement_covariance_estimator_outputs.end();
       ++i)
  {
    shifts.push_back((*i)->get_shift(inverse_incremental_temperature));
  }
  return shifts;
}

std::vector<arma::colvec> HMMEnsembleFactorVariables::get_deterministic_shifts() const
{
  std::vector<arma::colvec> shifts;
  shifts.reserve(this->measurement_covariance_estimator_outputs.size());
  for (auto i=this->measurement_covariance_estimator_outputs.begin();
       i!=this->measurement_covariance_estimator_outputs.end();
       ++i)
  {
    shifts.push_back((*i)->get_deterministic_shift());
  }
  return shifts;
}

/*
std::vector<arma::mat> HMMEnsembleFactorVariables::get_unconditional_measurement_covariances(const std::vector<arma::mat> &Cyys,
                                                                    double inverse_incremental_temperature) const
{
  std::vector<arma::mat> unconditional_measurement_covariances;
  unconditional_measurement_covariances.reserve(this->measurement_covariance_estimator_outputs.size());
  for (size_t i=0;
       i<this->measurement_covariance_estimator_outputs.size();
       ++i)
  {
    unconditional_measurement_covariances.push_back(this->measurement_covariance_estimator_outputs[i]->get_unconditional_measurement_covariance(Cyys[i],
                                                                                                                                                inverse_incremental_temperature));
  }
  return unconditional_measurement_covariances;
}
*/

std::vector<arma::colvec*> HMMEnsembleFactorVariables::get_measurements() const
{
  std::vector<arma::colvec*> measurements;
  measurements.reserve(this->measurement_covariance_estimator_outputs.size());
  for (auto i=this->measurement_covariance_estimator_outputs.begin();
       i!=this->measurement_covariance_estimator_outputs.end();
       ++i)
  {
    measurements.push_back((*i)->get_measurement());
  }
  return measurements;
}

double HMMEnsembleFactorVariables::evaluate_ensemble_likelihood_ratios(const Index* index,
                                                                       double inverse_incremental_temperature)
{
  double result = 0.0;
  for (size_t i=0;
       i<this->measurement_covariance_estimator_outputs.size();
       ++i)
  {
    result = result + this->measurement_covariance_estimator_outputs[i]->evaluate_ensemble_likelihood_ratio(inverse_incremental_temperature);
  }
  return result;
}

/*
double HMMEnsembleFactorVariables::evaluate_ensemble_likelihood_ratios(const Index* index,
                                                                       double inverse_incremental_temperature,
                                                                       const Parameters &conditioned_on_parameters)
{
  double result = 0.0;
  for (size_t i=0;
       i<this->measurement_covariance_estimator_outputs.size();
       ++i)
  {
    result = result + this->measurement_covariance_estimator_outputs[i]->evaluate_ensemble_likelihood_ratio(inverse_incremental_temperature,
                                                                                                            conditioned_on_parameters);
  }
  return result;
}
*/

double HMMEnsembleFactorVariables::subsample_evaluate_ensemble_likelihood_ratios(const Index* index,
                                                                                 double inverse_incremental_temperature)
{
  double result = 0.0;
  for (size_t i=0;
       i<this->measurement_covariance_estimator_outputs.size();
       ++i)
  {
    result = result + this->measurement_covariance_estimator_outputs[i]->subsample_evaluate_ensemble_likelihood_ratio(inverse_incremental_temperature);
  }
  return result;
}

/*
double HMMEnsembleFactorVariables::subsample_evaluate_ensemble_likelihood_ratios(const Index* index,
                                                                                 double inverse_incremental_temperature,
                                                                                 const Parameters &conditioned_on_parameters)
{
  double result = 0.0;
  for (size_t i=0;
       i<this->measurement_covariance_estimator_outputs.size();
       ++i)
  {
    result = result + this->measurement_covariance_estimator_outputs[i]->subsample_evaluate_ensemble_likelihood_ratio(inverse_incremental_temperature,
                                                                                                                      conditioned_on_parameters);
  }
  return result;
}
*/

double HMMEnsembleFactorVariables::evaluate_likelihoods(const Index* index)
{
  double result = 0.0;
  for (size_t i=0;
       i<this->measurement_covariance_estimator_outputs.size();
       ++i)
  {
    result = result + this->measurement_covariance_estimator_outputs[i]->evaluate_likelihood();
  }
  return result;
}

/*
double HMMEnsembleFactorVariables::evaluate_likelihoods(const Index* index,
                                                        const Parameters &conditioned_on_parameters)
{
  double result = 0.0;
  for (size_t i=0;
       i<this->measurement_covariance_estimator_outputs.size();
       ++i)
  {
    result = result + this->measurement_covariance_estimator_outputs[i]->evaluate_likelihood(conditioned_on_parameters);
  }
  return result;
}
*/

double HMMEnsembleFactorVariables::subsample_evaluate_likelihoods(const Index* index)
{
  double result = 0.0;
  for (size_t i=0;
       i<this->measurement_covariance_estimator_outputs.size();
       ++i)
  {
    result = result + this->measurement_covariance_estimator_outputs[i]->subsample_evaluate_likelihood();
  }
  return result;
}

/*
double HMMEnsembleFactorVariables::subsample_evaluate_likelihoods(const Index* index,
                                                                  const Parameters &conditioned_on_parameters)
{
  double result = 0.0;
  for (size_t i=0;
       i<this->measurement_covariance_estimator_outputs.size();
       ++i)
  {
    result = result + this->measurement_covariance_estimator_outputs[i]->subsample_evaluate_likelihood(conditioned_on_parameters);
  }
  return result;
}
*/

const EnsembleFactors* HMMEnsembleFactorVariables::get_ensemble_factors() const
{
  return this->hmm_ensemble_factors;
}

void HMMEnsembleFactorVariables::forget_you_were_already_written_to_file()
{
}

void HMMEnsembleFactorVariables::write_to_file(const std::string &directory_name,
                                               const std::string &index) const
{
  for (size_t i=0;
       i<this->measurement_covariance_estimator_outputs.size();
       ++i)
  {
    std::string factor_directory_name = directory_name + "/ensemble_factor" + toString(i);
    this->measurement_covariance_estimator_outputs[i]->write_to_file(factor_directory_name,
                                                                     index);
  }
}

void HMMEnsembleFactorVariables::close_ofstreams()
{
  for (size_t i=0;
       i<this->measurement_covariance_estimator_outputs.size();
       ++i)
  {
    this->measurement_covariance_estimator_outputs[i]->close_ofstreams();
  }
}

/*
void HMMEnsembleFactorVariables::evaluate_smcfixed_part_of_measurement_covariances(const Index* index)
{
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    if (this->measurement_covariance_estimator_outputs[*i]!=NULL)
    {
      this->measurement_covariance_estimator_outputs[*i]->evaluate_smcfixed_part(this->particle->parameters);
    }
  }
}

void HMMEnsembleFactorVariables::evaluate_smcfixed_part_of_measurement_covariances(const Index* index,
                                                                  const Parameters &conditioned_on_parameters)
{
  Parameters all_parameters = this->particle->parameters.merge(conditioned_on_parameters);
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    if (this->measurement_covariance_estimator_outputs[*i]!=NULL)
    {
      this->measurement_covariance_estimator_outputs[*i]->evaluate_smcfixed_part(all_parameters);
    }
  }
}

void HMMEnsembleFactorVariables::subsample_evaluate_smcfixed_part_of_measurement_covariances(const Index* index,
                                                                            const Parameters &conditioned_on_parameters)
{
  Parameters all_parameters = this->particle->parameters.merge(conditioned_on_parameters);
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    if (this->measurement_covariance_estimator_outputs[*i]!=NULL)
    {
      this->measurement_covariance_estimator_outputs[*i]->subsample_evaluate_smcfixed_part(all_parameters);
    }
  }
}

double HMMEnsembleFactorVariables::evaluate_smcadaptive_part_given_smcfixed_measurement_covariances(const Index* index)
{
  double result = 0.0;
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    if (this->measurement_covariance_estimator_outputs[*i]!=NULL)
    {
      this->measurement_covariance_estimator_outputs[*i]->evaluate_smcadaptive_part_given_smcfixed(this->particle->parameters);
      result = result + this->measurement_covariance_estimator_outputs[*i]->log_measurement_covariance;
    }
  }
  //this->target_evaluated = result;
  return result;
}

double HMMEnsembleFactorVariables::evaluate_smcadaptive_part_given_smcfixed_measurement_covariances(const Index* index,
                                                                                   const Parameters &conditioned_on_parameters)
{
  Parameters all_parameters = this->particle->parameters.merge(conditioned_on_parameters);
  double result = 0.0;
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    if (this->measurement_covariance_estimator_outputs[*i]!=NULL)
    {
      this->measurement_covariance_estimator_outputs[*i]->evaluate_smcadaptive_part_given_smcfixed(all_parameters);
      result = result + this->measurement_covariance_estimator_outputs[*i]->log_measurement_covariance;
    }
  }
  //this->target_evaluated = result;
  return result;
}

double HMMEnsembleFactorVariables::subsample_evaluate_smcadaptive_part_given_smcfixed_measurement_covariances(const Index* index,
                                                                                             const Parameters &conditioned_on_parameters)
{
  Parameters all_parameters = this->particle->parameters.merge(conditioned_on_parameters);
  double result = 0.0;
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    if (this->measurement_covariance_estimator_outputs[*i]!=NULL)
    {
      this->measurement_covariance_estimator_outputs[*i]->subsample_evaluate_smcadaptive_part_given_smcfixed(all_parameters);
      result = result + this->measurement_covariance_estimator_outputs[*i]->subsample_log_measurement_covariance;
    }
  }
  //this->subsample_target_evaluated = result;
  return result;
}

double HMMEnsembleFactorVariables::evaluate_measurement_covariances(const Index* index)
{
  double result = 0.0;
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    if (this->measurement_covariance_estimator_outputs[*i]!=NULL)
    {
      result = result + this->measurement_covariance_estimator_outputs[*i]->evaluate(this->particle->parameters);
    }
  }
  //this->target_evaluated = result;
  return result;
}

double HMMEnsembleFactorVariables::evaluate_measurement_covariances(const Index* index,
                                                   const Parameters &conditioned_on_parameters)
{
  Parameters all_parameters = this->particle->parameters.merge(conditioned_on_parameters);
  double result = 0.0;
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    if (this->measurement_covariance_estimator_outputs[*i]!=NULL)
    {
      result = result + this->measurement_covariance_estimator_outputs[*i]->evaluate(all_parameters);
    }
  }
  //this->target_evaluated = result;
  return result;
}

double HMMEnsembleFactorVariables::subsample_evaluate_measurement_covariances(const Index* index)
{
  double result = 0.0;
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    if (this->measurement_covariance_estimator_outputs[*i]!=NULL)
    {
      result = result + this->measurement_covariance_estimator_outputs[*i]->subsample_evaluate(this->particle->parameters);
    }
  }
  //this->subsample_target_evaluated = result;
  return result;
}

double HMMEnsembleFactorVariables::subsample_evaluate_measurement_covariances(const Index* index,
                                                             const Parameters &conditioned_on_parameters)
{
  Parameters all_parameters = this->particle->parameters.merge(conditioned_on_parameters);
  double result = 0.0;
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    if (this->measurement_covariance_estimator_outputs[*i]!=NULL)
    {
      result = result + this->measurement_covariance_estimator_outputs[*i]->subsample_evaluate(all_parameters);
    }
  }
  //this->subsample_target_evaluated = result;
  return result;
}

arma::mat HMMEnsembleFactorVariables::direct_get_gradient_of_log(const std::string &variable)
{
  arma::mat current_parameter = this->particle->parameters[variable];
  arma::mat result(current_parameter.n_rows,current_parameter.n_cols);
  for (std::vector<MeasurementCovarianceEstimatorOutput*>::const_iterator i=this->measurement_covariance_estimator_outputs.begin();
       i!=this->measurement_covariance_estimator_outputs.end();
       ++i)
  {
    if (*i!=NULL)
    {
      result = result + (*i)->get_gradient_of_log(variable,
                                                  this->particle->parameters);
    }
  }
  //this->target_gradients_of_log[variable] = result;
  return result;
}

arma::mat HMMEnsembleFactorVariables::direct_get_gradient_of_log(const std::string &variable,
                                               const Parameters &conditioned_on_parameters)
{
  Parameters all_parameters = this->particle->parameters.merge(conditioned_on_parameters);
  arma::mat current_parameter = this->particle->parameters[variable];
  arma::mat result(current_parameter.n_rows,current_parameter.n_cols);
  for (std::vector<MeasurementCovarianceEstimatorOutput*>::const_iterator i=this->measurement_covariance_estimator_outputs.begin();
       i!=this->measurement_covariance_estimator_outputs.end();
       ++i)
  {
    if (*i!=NULL)
    {
      result = result + (*i)->get_gradient_of_log(variable,
                                                  all_parameters);
    }
  }
  //this->target_gradients_of_log[variable] = result;
  return result;
}

arma::mat HMMEnsembleFactorVariables::direct_subsample_get_gradient_of_log(const std::string &variable)
{
  arma::mat current_parameter = this->particle->parameters[variable];
  arma::mat result(current_parameter.n_rows,current_parameter.n_cols);
  for (std::vector<MeasurementCovarianceEstimatorOutput*>::const_iterator i=this->measurement_covariance_estimator_outputs.begin();
       i!=this->measurement_covariance_estimator_outputs.end();
       ++i)
  {
    if (*i!=NULL)
    {
      result = result + (*i)->subsample_get_gradient_of_log(variable,
                                                            this->particle->parameters);
    }
  }
  //this->subsample_target_gradients_of_log[variable] = result;
  return result;
}

arma::mat HMMEnsembleFactorVariables::direct_subsample_get_gradient_of_log(const std::string &variable,
                                                                      const Parameters &conditioned_on_parameters)
{
  Parameters all_parameters = this->particle->parameters.merge(conditioned_on_parameters);
  arma::mat current_parameter = this->particle->parameters[variable];
  arma::mat result(current_parameter.n_rows,current_parameter.n_cols);
  for (std::vector<MeasurementCovarianceEstimatorOutput*>::const_iterator i=this->measurement_covariance_estimator_outputs.begin();
       i!=this->measurement_covariance_estimator_outputs.end();
       ++i)
  {
    if (*i!=NULL)
    {
      result = result + (*i)->subsample_get_gradient_of_log(variable,
                                                            all_parameters);
    }
  }
  //this->subsample_target_gradients_of_log[variable] = result;
  return result;
}

arma::mat HMMEnsembleFactorVariables::direct_get_gradient_of_log(const Index* index,
                                                            const std::string &variable)
{
  arma::mat current_parameter = this->particle->parameters[variable];
  arma::mat result(current_parameter.n_rows,current_parameter.n_cols);
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    if (this->measurement_covariance_estimator_outputs[*i]!=NULL)
    {
      result = result + this->measurement_covariance_estimator_outputs[*i]->get_gradient_of_log(variable,
                                                                                    this->particle->parameters);
    }
  }
  //this->target_gradients_of_log[variable] = result;
  return result;
}

arma::mat HMMEnsembleFactorVariables::direct_get_gradient_of_log(const Index* index,
                                                            const std::string &variable,
                                                            const Parameters &conditioned_on_parameters)
{
  Parameters all_parameters = this->particle->parameters.merge(conditioned_on_parameters);
  arma::mat current_parameter = this->particle->parameters[variable];
  arma::mat result(current_parameter.n_rows,current_parameter.n_cols);
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    if (this->measurement_covariance_estimator_outputs[*i]!=NULL)
    {
      result = result + this->measurement_covariance_estimator_outputs[*i]->get_gradient_of_log(variable,
                                                                                    all_parameters);
    }
  }
  //this->target_gradients_of_log[variable] = result;
  return result;
}

arma::mat HMMEnsembleFactorVariables::direct_subsample_get_gradient_of_log(const Index* index,
                                                                      const std::string &variable)
{
  arma::mat current_parameter = this->particle->parameters[variable];
  arma::mat result(current_parameter.n_rows,current_parameter.n_cols);
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    if (this->measurement_covariance_estimator_outputs[*i]!=NULL)
    {
      result = result + this->measurement_covariance_estimator_outputs[*i]->subsample_get_gradient_of_log(variable,
                                                                                              this->particle->parameters);
    }
  }
  //this->subsample_target_gradients_of_log[variable] = result;
  return result;
}

arma::mat HMMEnsembleFactorVariables::direct_subsample_get_gradient_of_log(const Index* index,
                                                                      const std::string &variable,
                                                                      const Parameters &conditioned_on_parameters)
{
  Parameters all_parameters = this->particle->parameters.merge(conditioned_on_parameters);
  arma::mat current_parameter = this->particle->parameters[variable];
  arma::mat result(current_parameter.n_rows,current_parameter.n_cols);
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    if (this->measurement_covariance_estimator_outputs[*i]!=NULL)
    {
      result = result + this->measurement_covariance_estimator_outputs[*i]->subsample_get_gradient_of_log(variable,
                                                                                              all_parameters);
    }
  }
  //this->subsample_target_gradients_of_log[variable] = result;
  return result;
}
*/
