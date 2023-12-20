// Include the class definition.
#include "ensemble.h"
#include "ensemble_factor_variables.h"
#include "ensemble_factors.h"
#include "ensemble_kalman.h"
#include "move_output.h"
#include "single_point_move_output.h"
#include "utils.h"

// Include other class definitions that are needed.

// Public functions //

// Default constructor and destructor.

// Comment about function.
Ensemble::Ensemble()
{
	// Set all the pointers owned by this class to be something (possibly NULL).
  //this->packing_instructions = NULL;
  this->ensemble_factors = NULL;
  this->log_normalising_constant_ratio = 0.0;
}

// Comment about function.
Ensemble::~Ensemble()
{
	// Delete all the pointers owned by this class.
  for (std::vector<MoveOutput*>::iterator i=this->members.begin();
       i!=this->members.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  
  // Delete all the pointers owned by this class.
  for (std::vector<MoveOutput*>::iterator i=this->predicted_members.begin();
       i!=this->predicted_members.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
}

// Other constuctors.

Ensemble::Ensemble(EnsembleFactors* ensemble_factors_in)
{
  this->ensemble_factors = ensemble_factors_in;
}

Ensemble::Ensemble(std::vector<Parameters> &initial_values_in,
                   EnsembleFactors* factors_in)
{
  this->setup(initial_values_in,factors_in);
}

void Ensemble::setup(std::vector<Parameters> &initial_values_in,
                     EnsembleFactors* factors_in)
{
  size_t number_of_particles_in = initial_values_in.size();
  this->members.reserve(number_of_particles_in);
  this->unnormalised_log_weights = arma::colvec(number_of_particles_in);
  this->log_normalising_constant_ratio = 0.0;
  
  size_t counter = 0;
  for (std::vector<Parameters>::iterator i=initial_values_in.begin();
       i!=initial_values_in.end();
       ++i, ++counter)
  {
    this->push_back(std::move(*i),factors_in);
  }
  
  this->ensemble_factors = factors_in;
}

// Everything you need to copy the class.

// The copy constructor.
Ensemble::Ensemble(const Ensemble &another)
{
	this->make_copy(another);
}

Ensemble::Ensemble(Ensemble &&another)
{
  this->make_copy(std::move(another));
}

// Returns a pointer to a base class of type Ensemble.
Ensemble* Ensemble::duplicate() const
{
	return new Ensemble(*this);
}

// Copy all the members of the class.
void Ensemble::make_copy(const Ensemble &another)
{
	// Copy all members, duplicating the memory where appropriate.
  
  this->members.resize(0);
  this->members.reserve(another.members.size());
  for (std::vector<MoveOutput*>::const_iterator i=another.members.begin();
       i!=another.members.end();
       ++i)
  {
    if (*i!=NULL)
      this->members.push_back((*i)->duplicate());
    else
      this->members.push_back(NULL);
  }
  
  this->predicted_members.resize(0);
  this->predicted_members.reserve(another.predicted_members.size());
  for (std::vector<MoveOutput*>::const_iterator i=another.predicted_members.begin();
       i!=another.predicted_members.end();
       ++i)
  {
    if (*i!=NULL)
      this->predicted_members.push_back((*i)->duplicate());
    else
      this->predicted_members.push_back(NULL);
  }
  
  //this->members = another.members;
  this->packed_members = another.packed_members;
  //this->predicted_members = another.predicted_members;
  this->partially_packed_members_row = another.partially_packed_members_row;
  this->partially_packed_members_col = another.partially_packed_members_col;
  this->partially_packed_predicted_members_col = another.partially_packed_predicted_members_col;
  //this->packing_instructions = another.packing_instructions;
  this->packed_measurement_states = another.packed_measurement_states;
  this->partially_packed_measurement_states = another.partially_packed_measurement_states;
  this->Cxys = another.Cxys;
  this->Cyys = another.Cyys;
  this->myys = another.myys;
  this->kalman_gains = another.kalman_gains;
  this->inv_Cxx = another.inv_Cxx;
  //this->measurements = another.measurements;
  this->ensemble_factors = another.ensemble_factors;
  this->log_normalising_constant_ratio = another.log_normalising_constant_ratio;
  this->unnormalised_log_weights = another.unnormalised_log_weights;
  this->schedule_parameters = another.schedule_parameters;
  this->ess = another.ess;
  //this->partially_packed_measurement_random_shifts = another.partially_packed_measurement_random_shifts;
}

void Ensemble::make_copy(Ensemble &&another)
{
  // Copy all members, duplicating the memory where appropriate.
  
  this->members.resize(0);
  this->members.reserve(another.members.size());
  for (std::vector<MoveOutput*>::const_iterator i=another.members.begin();
       i!=another.members.end();
       ++i)
  {
    if (*i!=NULL)
      this->members.push_back(*i);
    else
      this->members.push_back(NULL);
  }
  
  this->predicted_members.resize(0);
  this->predicted_members.reserve(another.predicted_members.size());
  for (std::vector<MoveOutput*>::const_iterator i=another.predicted_members.begin();
       i!=another.predicted_members.end();
       ++i)
  {
    if (*i!=NULL)
      this->predicted_members.push_back(*i);
    else
      this->predicted_members.push_back(NULL);
  }
  
  //this->members = another.members;
  this->packed_members = std::move(another.packed_members);
  //this->predicted_members = std::move(another.predicted_members);
  this->partially_packed_members_row = std::move(another.partially_packed_members_row);
  this->partially_packed_members_col = std::move(another.partially_packed_members_col);
  this->partially_packed_predicted_members_col = std::move(another.partially_packed_predicted_members_col);
  //this->packing_instructions = another.packing_instructions;
  this->packed_measurement_states = std::move(another.packed_measurement_states);
  this->partially_packed_measurement_states = std::move(another.partially_packed_measurement_states);
  this->Cxys = std::move(another.Cxys);
  this->Cyys = std::move(another.Cyys);
  this->myys = std::move(another.myys);
  this->kalman_gains = std::move(another.kalman_gains);
  this->inv_Cxx = std::move(another.inv_Cxx);
  //this->measurements = another.measurements;
  this->ensemble_factors = another.ensemble_factors;
  this->log_normalising_constant_ratio = std::move(another.log_normalising_constant_ratio);
  this->unnormalised_log_weights = std::move(another.unnormalised_log_weights);
  this->schedule_parameters = std::move(another.schedule_parameters);
  this->ess = std::move(another.ess);
  //this->partially_packed_measurement_random_shifts = another.partially_packed_measurement_random_shifts;
  
  another.members = std::vector<MoveOutput*>();

  //another.members = another.members;
  another.packed_members = arma::mat();
  another.predicted_members = std::vector<MoveOutput*>();
  another.partially_packed_members_row = std::vector<arma::rowvec>();
  another.partially_packed_members_col = std::vector<arma::colvec>();
  another.partially_packed_predicted_members_col = std::vector<arma::colvec>();
  //another.packing_instructions = another.packing_instructions;
  another.packed_measurement_states = std::vector<arma::mat>();
  another.partially_packed_measurement_states = std::vector< std::vector<arma::rowvec> >();
  another.Cxys = std::vector<arma::mat>();
  another.Cyys = std::vector<arma::mat>();
  another.myys = std::vector<arma::colvec>();
  another.kalman_gains = std::vector<arma::mat>();
  another.inv_Cxx = arma::mat();
  //another.measurements = another.measurements;
  another.ensemble_factors = NULL;
  another.log_normalising_constant_ratio = 0.0;
  another.unnormalised_log_weights = arma::colvec();
  another.schedule_parameters = Parameters();
  another.ess = 0.0;
}

// The = operator.
Ensemble& Ensemble::operator=(const Ensemble &another)
{
	if(this==&another)//a==a
		return *this;
  
  for (std::vector<MoveOutput*>::iterator i=this->members.begin();
       i!=this->members.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  
  // Delete all the pointers owned by this class.
  for (std::vector<MoveOutput*>::iterator i=this->predicted_members.begin();
       i!=this->predicted_members.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  
  this->members.clear();
  this->predicted_members.clear();
	
	this->make_copy(another);
  
  return *this;
}

Ensemble& Ensemble::operator=(Ensemble &&another)
{
  if(this==&another)//a==a
    return *this;
  
  for (std::vector<MoveOutput*>::iterator i=this->members.begin();
       i!=this->members.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  
  // Delete all the pointers owned by this class.
  for (std::vector<MoveOutput*>::iterator i=this->predicted_members.begin();
       i!=this->predicted_members.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  
  this->members.clear();
  this->predicted_members.clear();
  
  this->make_copy(std::move(another));
  
  return *this;
}

void Ensemble::reserve(size_t number_of_ensemble_members_in)
{
  this->members.reserve(number_of_ensemble_members_in);
  this->unnormalised_log_weights = arma::colvec(number_of_ensemble_members_in);
}

void Ensemble::push_back(Particle &&particle_in)
{
  MoveOutput* single_particle = new SinglePointMoveOutput(std::move(particle_in));
  this->members.push_back(single_particle);
}

void Ensemble::push_back(Parameters &&parameters_in,
                         EnsembleFactors* factors_in)
{
  MoveOutput* single_particle = new SinglePointMoveOutput(std::move(parameters_in),
                                                          factors_in);
  this->members.push_back(single_particle);
}

void Ensemble::push_back(MoveOutput* move_output_in)
{
  this->members.push_back(move_output_in);
}

Particle* Ensemble::add_ensemble_member()
{
  MoveOutput* single_particle = new SinglePointMoveOutput();
  this->members.push_back(single_particle);
  return &this->members.back()->back();
}

MoveOutput* Ensemble::operator[](const size_t &i)
{
  return(this->members[i]);
}

MoveOutput* Ensemble::operator[](const size_t &i) const
{
  return(this->members[i]);
}

MoveOutput* Ensemble::back()
{
  return this->members.back();
}

MoveOutput* Ensemble::back() const
{
  return this->members.back();
}

size_t Ensemble::size() const
{
  return(this->members.size());
}

void Ensemble::update_weights(const arma::colvec &latest_unnormalised_log_incremental_weights)
{
  this->unnormalised_log_weights = latest_unnormalised_log_incremental_weights;
}

double Ensemble::calculate_log_normalising_constant()
{
  this->log_normalising_constant_ratio = this->ensemble_factors->get_incremental_likelihood(this);
  this->ess = exp(2.0*log_sum_exp(this->unnormalised_log_weights) - log_sum_exp(2.0*this->unnormalised_log_weights));
  return this->log_normalising_constant_ratio;
}

double Ensemble::calculate_inversion_log_normalising_constant(double inverse_incremental_temperature)
{
  this->log_normalising_constant_ratio = this->ensemble_factors->get_inversion_incremental_likelihood(this,
                                                                                                      inverse_incremental_temperature);
  this->ess = exp(2.0*log_sum_exp(this->unnormalised_log_weights) - log_sum_exp(2.0*this->unnormalised_log_weights));
  return this->log_normalising_constant_ratio;
}

arma::rowvec Ensemble::get_output_lengths() const
{
  arma::rowvec output_lengths(this->size());
  
  for (size_t i=0; i<this->members.size(); ++i)
  {
    output_lengths[i] = this->members[i]->length();
  }
  
  return output_lengths;
}

void Ensemble::find_measurement_covariances()
{
  this->Cxys.clear();
  this->myys.clear();
  this->Cyys.clear();
  this->Cxys.reserve(this->packed_measurement_states.size());
  this->Cyys.reserve(this->packed_measurement_states.size());
  
  for (size_t i=0;
       i<this->packed_measurement_states.size();
       ++i)
  {
    this->Cxys.push_back(arma::cov(this->packed_members,this->packed_measurement_states[i]));
    //this->Cyys.push_back(arma::cov(this->packed_measurement_states[i] + inverse_incremental_temperature*measurement_covariances[i]));
    this->Cyys.push_back(arma::cov(this->packed_measurement_states[i]));
    this->myys.push_back(arma::conv_to<arma::colvec>::from(arma::mean(this->packed_measurement_states[i],0)));
  }
  
  if (this->ensemble_factors->need_Cxx()==true)
  {
    this->inv_Cxx = arma::inv_sympd(arma::cov(this->packed_members));
    this->ensemble_factors->find_Cygivenx(this->inv_Cxx,
                                          this->Cxys,
                                          this->Cyys);
  }
  //this->ensemble_factors->find_reduced_Cygivenx(this->Cxys,
  //                                              this->Cyys,
  //                                              this->packed_members);
  
}

void Ensemble::find_measurement_covariances(const Parameters &conditioned_on_parameters)
{
  this->Cxys.clear();
  this->myys.clear();
  this->Cyys.clear();
  this->Cxys.reserve(this->packed_measurement_states.size());
  this->Cyys.reserve(this->packed_measurement_states.size());
  
  for (size_t i=0;
       i<this->packed_measurement_states.size();
       ++i)
  {
    this->Cxys.push_back(arma::cov(this->packed_members,this->packed_measurement_states[i]));
    //this->Cyys.push_back(arma::cov(this->packed_measurement_states[i] + measurement_covariances[i]));
    this->Cyys.push_back(arma::cov(this->packed_measurement_states[i]));
    this->myys.push_back(arma::conv_to<arma::colvec>::from(arma::mean(this->packed_measurement_states[i],0)));
  }
  
  if (this->ensemble_factors->need_Cxx()==true)
  {
    this->inv_Cxx = arma::inv_sympd(arma::cov(this->packed_members));
    this->ensemble_factors->find_Cygivenx(this->inv_Cxx,
                                          this->Cxys,
                                          this->Cyys);
  }
}

/*
void Ensemble::pack_measurements()
{
  this->measurements = this->ensemble_factors->pack_measurements();
}
*/

void Ensemble::set_temperature(double temperature_in)
{
  Rcpp::stop("Ensemble::set_temperature - should no longer be used.");
  //this->ensemble_factors->set_temperature(temperature_in);
}

double Ensemble::get_inverse_incremental_temperature() const
{
  Rcpp::stop("Ensemble::get_inverse_incremental_temperature - should no longer be used.");
  //return this->ensemble_factors->get_inverse_incremental_temperature();
}

void Ensemble::precompute_gaussian_covariance(double inverse_incremental_temperature)
{
  ensemble_factors->precompute_gaussian_covariance(inverse_incremental_temperature);
}

arma::mat Ensemble::get_packed_members() const
{
  return this->packed_members;
}

void Ensemble::close_ofstreams()
{
  for (size_t i=0; i<this->members.size(); ++i)
  {
    this->members[i]->close_ofstreams();
  }
}

/*
void Ensemble::set_previous_target_evaluated_to_target_evaluated()
{
  for (std::vector<MoveOutput*>::iterator i = this->members.begin();
       i != this->members.end();
       ++i)
  {
    (*i)->back().previous_target_evaluated = (*i)->back().target_evaluated;
    (*i)->back().previous_ensemble_target_evaluated = (*i)->back().ensemble_target_evaluated;
    //(*i)->back().previous_target_gradients_of_log = (*i)->back().target_gradients_of_log;
  }
}

void Ensemble::subsample_set_previous_target_evaluated_to_target_evaluated()
{
  for (std::vector<MoveOutput*>::iterator i = this->members.begin();
       i != this->members.end();
       ++i)
  {
    (*i)->back().subsample_previous_target_evaluated = (*i)->back().subsample_target_evaluated;
    (*i)->back().subsample_previous_ensemble_target_evaluated = (*i)->back().subsample_ensemble_target_evaluated;
    //(*i)->back().subsample_previous_target_gradients_of_log = (*i)->back().subsample_target_gradients_of_log;
  }
}
*/
