// Include the class definition.
#include "ensemble.h"
#include "ensemble_factor_variables.h"
#include "ensemble_factors.h"
#include "ensemble_kalman.h"
#include "move_output.h"
#include "single_point_move_output.h"

// Include other class definitions that are needed.

// Public functions //

// Default constructor and destructor.

// Comment about function.
Ensemble::Ensemble()
{
	// Set all the pointers owned by this class to be something (possibly NULL).
  this->packing_instructions = NULL;
  this->ensemble_factors = NULL;
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
}

// Other constuctors.

// Everything you need to copy the class.

// The copy constructor.
Ensemble::Ensemble(const Ensemble &another)
{
	this->make_copy(another);
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
  
  //this->members = another.members;
  this->packed_members = another.packed_members;
  this->partially_packed_members_row = another.partially_packed_members_row;
  this->partially_packed_members_col = another.partially_packed_members_col;
  this->packing_instructions = another.packing_instructions;
  this->packed_measurement_states = another.packed_measurement_states;
  this->partially_packed_measurement_states = another.partially_packed_measurement_states;
  this->Cxys = another.Cxys;
  this->Cyys = another.Cyys;
  this->inv_Cxx = another.inv_Cxx;
  this->ensemble_factors = another.ensemble_factors;
  //this->partially_packed_measurement_random_shifts = another.partially_packed_measurement_random_shifts;
}

// The = operator.
void Ensemble::operator=(const Ensemble &another)
{
	if(this==&another)//a==a
		return;
  
  for (std::vector<MoveOutput*>::iterator i=this->members.begin();
       i!=this->members.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  
  this->members.clear();
	
	this->make_copy(another);
}

void Ensemble::reserve(size_t number_of_ensemble_members_in)
{
  this->members.reserve(number_of_ensemble_members_in);
}

void Ensemble::push_back(const Particle &particle_in)
{
  MoveOutput* single_particle = new SinglePointMoveOutput(particle_in);
  this->members.push_back(single_particle);
}

void Ensemble::push_back(MoveOutput* move_output_in)
{
  this->members.push_back(move_output_in);
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

void Ensemble::find_measurement_covariances()
{
  this->Cxys.clear();
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
  }
  
}

void Ensemble::find_measurement_covariances(const Parameters &conditioned_on_parameters)
{
  this->Cxys.clear();
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
  }
  
  if (this->ensemble_factors->need_Cxx()==true)
  {
    this->inv_Cxx = arma::cov(this->packed_members);
    this->ensemble_factors->find_Cygivenx(this->inv_Cxx,
                                          this->Cxys,
                                          this->Cyys);
  }
}

void Ensemble::set_temperature(double temperature_in)
{
  this->ensemble_factors->set_temperature(temperature_in);
}

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
