#include "particles.h"
#include "sequential_smc_worker.h"
#include "smc.h"
#include "utils.h"
#include "move_output.h"
#include "single_point_move_output.h"

//#include <array>

Particles::Particles()
{
  this->resampled_flag = false;
}

Particles::Particles(size_t number_of_particles_in)
{
  this->reserve(number_of_particles_in);
}

Particles::Particles(const std::vector< MoveOutput* > &particles_in)
{
  this->particles = particles_in;
  this->resampling_variables = arma::colvec(this->particles.size());
  this->ancestor_variables = std::vector<size_t>(this->particles.size());
  this->unnormalised_log_weights = arma::colvec(this->particles.size());
  this->normalised_log_weights = arma::colvec(this->particles.size());
  this->previous_normalised_log_weights = arma::colvec(this->particles.size());
  this->incremental_log_weights = arma::colvec(this->particles.size());
  this->log_normalising_constant_ratio = 0.0;
  this->resampled_flag = false;
}

Particles::Particles(std::vector<Parameters> &initial_values_in,
                     const arma::colvec &log_probabilities_of_initial_values_in,
                     Factors* factors_in,
                     const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
                     const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in)
{
  this->setup(initial_values_in,
              log_probabilities_of_initial_values_in,
              factors_in,
              proposals_to_transform_for_in,
              proposals_to_find_gradient_for_in);
}

void Particles::setup(std::vector<Parameters> &initial_values_in,
                      const arma::colvec &log_probabilities_of_initial_values_in,
                      Factors* factors_in,
                      const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
                      const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in)
{
  if (initial_values_in.size()!=log_probabilities_of_initial_values_in.n_rows)
  {
    Rcpp::stop("Particles(initial_values,probs) constructor: values and probabilities need to be of the same length.");
  }
  
  if (factors_in==NULL)
  {
    Rcpp::stop("Particles(initial_values,probs) constructor: factors need to be non-null.");
  }
  
  size_t number_of_particles_in = initial_values_in.size();
  this->particles.reserve(number_of_particles_in);
  this->resampling_variables = arma::colvec(number_of_particles_in);
  this->ancestor_variables = std::vector<size_t>(number_of_particles_in);
  this->unnormalised_log_weights = arma::colvec(number_of_particles_in);
  this->normalised_log_weights = arma::colvec(number_of_particles_in);
  this->previous_normalised_log_weights = arma::colvec(number_of_particles_in);
  this->incremental_log_weights = arma::colvec(number_of_particles_in);
  this->log_normalising_constant_ratio = 0.0;
  this->resampled_flag = false;
  
  size_t counter = 0;
  for (std::vector<Parameters>::iterator i=initial_values_in.begin();
       i!=initial_values_in.end();
       ++i, ++counter)
  {
    //arma::mat tau = (*i)["tau"];
    this->push_back(std::move(*i),factors_in,proposals_to_transform_for_in,proposals_to_find_gradient_for_in);
    //this->push_back(*i,factors_in,proposals_to_transform_for_in,proposals_to_find_gradient_for_in);
    this->particles.back()->back().previous_target_evaluated = log_probabilities_of_initial_values_in[counter];
  }
}

Particles::Particles(std::vector<Parameters> &initial_values_in,
                     const arma::colvec &log_probabilities_of_initial_values_in,
                     Factors* factors_in,
                     const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
                     const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in,
                     const Parameters &conditioned_on_parameters)
{
  this->setup(initial_values_in,
              log_probabilities_of_initial_values_in,
              factors_in,
              proposals_to_transform_for_in,
              proposals_to_find_gradient_for_in,
              conditioned_on_parameters);
}

Particles::Particles(std::vector<Parameters> &initial_values_in,
                     const arma::colvec &log_probabilities_of_initial_values_in,
                     Factors* factors_in,
                     const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
                     const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in,
                     const Parameters &conditioned_on_parameters,
                     const Parameters &sequencer_parameters)
{
  this->setup(initial_values_in,
              log_probabilities_of_initial_values_in,
              factors_in,
              proposals_to_transform_for_in,
              proposals_to_find_gradient_for_in,
              conditioned_on_parameters,
              sequencer_parameters);
}

void Particles::setup(std::vector<Parameters> &initial_values_in,
                      const arma::colvec &log_probabilities_of_initial_values_in,
                      Factors* factors_in,
                      const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
                      const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in,
                      const Parameters &conditioned_on_parameters)
{
  if (initial_values_in.size()!=log_probabilities_of_initial_values_in.n_rows)
  {
    Rcpp::stop("Particles(initial_values,probs) constructor: values and probabilities need to be of the same length.");
  }
  
  if (factors_in==NULL)
  {
    Rcpp::stop("Particles(initial_values,probs) constructor: factors need to be non-null.");
  }
  
  size_t number_of_particles_in = initial_values_in.size();
  this->particles.reserve(number_of_particles_in);
  this->resampling_variables = arma::colvec(number_of_particles_in);
  this->ancestor_variables = std::vector<size_t>(number_of_particles_in);
  this->unnormalised_log_weights = arma::colvec(number_of_particles_in);
  this->normalised_log_weights = arma::colvec(number_of_particles_in);
  this->previous_normalised_log_weights = arma::colvec(number_of_particles_in);
  this->incremental_log_weights = arma::colvec(number_of_particles_in);
  this->log_normalising_constant_ratio = 0.0;
  this->resampled_flag = false;
  
  size_t counter = 0;
  for (std::vector<Parameters>::iterator i=initial_values_in.begin();
       i!=initial_values_in.end();
       ++i, ++counter)
  {
    //arma::mat x = (*i)["x"];
    //std::cout << x << std::endl;
    i->merge_with_fixed(conditioned_on_parameters);
    this->push_back(std::move(*i),factors_in,proposals_to_transform_for_in,proposals_to_find_gradient_for_in);
    
    //std::cout << this->back()->back().parameters["x"] << std::endl;
    //std::cout << this->back()->back().parameters << std::endl;
    //this->push_back(*i,factors_in,proposals_to_transform_for_in,proposals_to_find_gradient_for_in);
    this->particles.back()->back().previous_target_evaluated = log_probabilities_of_initial_values_in[counter];
  }
}

void Particles::setup(std::vector<Parameters> &initial_values_in,
                      const arma::colvec &log_probabilities_of_initial_values_in,
                      Factors* factors_in,
                      const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
                      const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in,
                      const Parameters &conditioned_on_parameters,
                      const Parameters &sequencer_parameters)
{
  if (initial_values_in.size()!=log_probabilities_of_initial_values_in.n_rows)
  {
    Rcpp::stop("Particles(initial_values,probs) constructor: values and probabilities need to be of the same length.");
  }
  
  if (factors_in==NULL)
  {
    Rcpp::stop("Particles(initial_values,probs) constructor: factors need to be non-null.");
  }
  
  size_t number_of_particles_in = initial_values_in.size();
  this->particles.reserve(number_of_particles_in);
  this->resampling_variables = arma::colvec(number_of_particles_in);
  this->ancestor_variables = std::vector<size_t>(number_of_particles_in);
  this->unnormalised_log_weights = arma::colvec(number_of_particles_in);
  this->normalised_log_weights = arma::colvec(number_of_particles_in);
  this->previous_normalised_log_weights = arma::colvec(number_of_particles_in);
  this->incremental_log_weights = arma::colvec(number_of_particles_in);
  this->log_normalising_constant_ratio = 0.0;
  this->resampled_flag = false;
  
  size_t counter = 0;
  for (std::vector<Parameters>::iterator i=initial_values_in.begin();
       i!=initial_values_in.end();
       ++i, ++counter)
  {
    //arma::mat tau = (*i)["tau"];
    i->merge_with_fixed(conditioned_on_parameters);
    i->merge_with_fixed(sequencer_parameters);
    this->push_back(std::move(*i),factors_in,proposals_to_transform_for_in,proposals_to_find_gradient_for_in);
    this->particles.back()->back().previous_target_evaluated = log_probabilities_of_initial_values_in[counter];
  }
}

Particles::Particles(const Particles &another)
{
  this->make_copy(another);
}

Particles::Particles(Particles &&another)
{
  this->make_copy(std::move(another));
}

Particles::~Particles()
{
  for (std::vector<MoveOutput*>::iterator i=this->particles.begin();
       i!=this->particles.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
}

Particles& Particles::operator=(const Particles &another)
{
  if(this == &another)
    return *this;
  
  for (std::vector<MoveOutput*>::iterator i=this->particles.begin();
       i!=this->particles.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  
  this->particles.clear();

  this->make_copy(another);
  
  return *this;
}

Particles& Particles::operator=(Particles &&another)
{
  if(this == &another)
    return *this;
  
  for (std::vector<MoveOutput*>::iterator i=this->particles.begin();
       i!=this->particles.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  
  this->particles.clear();
  
  this->make_copy(std::move(another));
  
  return *this;
}

void Particles::reserve(size_t number_of_particles_in)
{
  this->particles.reserve(number_of_particles_in);
  this->resampling_variables = arma::colvec(number_of_particles_in);
  this->ancestor_variables = std::vector<size_t>(number_of_particles_in);
  this->unnormalised_log_weights = arma::colvec(number_of_particles_in);
  this->normalised_log_weights = arma::colvec(number_of_particles_in);
  this->previous_normalised_log_weights = arma::colvec(number_of_particles_in);
  this->incremental_log_weights = arma::colvec(number_of_particles_in);
  this->log_normalising_constant_ratio = 0.0;
}

/*
void Particles::push_back(const Parameters &parameters_in,
                          Factors* factors_in)
{
  MoveOutput* single_particle = new SinglePointMoveOutput(parameters_in,
                                                          factors_in);
  this->particles.push_back(single_particle);
}

void Particles::push_back(const Particle &particle_in)
{
  MoveOutput* single_particle = new SinglePointMoveOutput(particle_in);
  this->particles.push_back(single_particle);
}
*/

void Particles::push_back(Parameters &&parameters_in,
                          Factors* factors_in,
                          const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
                          const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in)
{
  //arma::mat tau = parameters_in["tau"];
  MoveOutput* single_particle = new SinglePointMoveOutput(std::move(parameters_in),
                                                          factors_in,
                                                          proposals_to_transform_for_in,
                                                          proposals_to_find_gradient_for_in);
  this->particles.push_back(single_particle);
}

void Particles::push_back(Particle &&particle_in)
{
  MoveOutput* single_particle = new SinglePointMoveOutput(std::move(particle_in));
  this->particles.push_back(single_particle);
}

void Particles::push_back(const Parameters &parameters_in,
                          Factors* factors_in,
                          const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
                          const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in)
{
  //arma::mat tau = parameters_in["tau"];
  MoveOutput* single_particle = new SinglePointMoveOutput(parameters_in,
                                                          factors_in,
                                                          proposals_to_transform_for_in,
                                                          proposals_to_find_gradient_for_in);
  this->particles.push_back(single_particle);
}

void Particles::push_back(const Particle &particle_in)
{
  MoveOutput* single_particle = new SinglePointMoveOutput(particle_in);
  this->particles.push_back(single_particle);
}

void Particles::push_back(MoveOutput* move_output_in)
{
  this->particles.push_back(move_output_in);
}

Particle* Particles::add_particle()
{
  MoveOutput* single_particle = new SinglePointMoveOutput();
  this->particles.push_back(single_particle);
  return &this->particles.back()->back();
}

/*
void Particles::push_back(const std::deque<Particle> &particle_in)
{
  this->particles.push_back((particle_in));
}
 */

void Particles::simulate_resampling_variables(RandomNumberGenerator &rng)
{
  for (arma::colvec::iterator i = this->resampling_variables.begin();
       i!=this->resampling_variables.end();
       ++i)
  {
    (*i) = runif(rng);
  }
}

void Particles::make_copy(const Particles &another)
{
  this->particles.resize(0);
  this->particles.reserve(another.particles.size());
  for (std::vector<MoveOutput*>::const_iterator i=another.particles.begin();
       i!=another.particles.end();
       ++i)
  {
    if (*i!=NULL)
      this->particles.push_back((*i)->duplicate());
    else
      this->particles.push_back(NULL);
  }
  
  this->resampling_variables = another.resampling_variables;
  this->ancestor_variables = another.ancestor_variables;
  this->unnormalised_log_weights = another.unnormalised_log_weights;
  this->normalised_log_weights = another.normalised_log_weights;
  this->previous_normalised_log_weights = another.previous_normalised_log_weights;
  this->incremental_log_weights = another.incremental_log_weights;
  this->log_normalising_constant_ratio = another.log_normalising_constant_ratio;
  this->resampled_flag = another.resampled_flag;
  this->schedule_parameters = another.schedule_parameters;
  this->ess = another.ess;
}

void Particles::make_copy(Particles &&another)
{
  this->particles.resize(0);
  this->particles.reserve(another.particles.size());
  for (std::vector<MoveOutput*>::const_iterator i=another.particles.begin();
       i!=another.particles.end();
       ++i)
  {
    if (*i!=NULL)
      this->particles.push_back(*i);
    else
      this->particles.push_back(NULL);
  }
  
  this->resampling_variables = std::move(another.resampling_variables);
  this->ancestor_variables = std::move(another.ancestor_variables);
  this->unnormalised_log_weights = std::move(another.unnormalised_log_weights);
  this->normalised_log_weights = std::move(another.normalised_log_weights);
  this->previous_normalised_log_weights = std::move(another.previous_normalised_log_weights);
  this->incremental_log_weights = std::move(another.incremental_log_weights);
  this->log_normalising_constant_ratio = std::move(another.log_normalising_constant_ratio);
  this->resampled_flag = std::move(another.resampled_flag);
  this->schedule_parameters = std::move(another.schedule_parameters);
  this->ess = std::move(another.ess);
  
  another.particles = std::vector<MoveOutput*>();
  another.resampling_variables = arma::colvec();
  another.ancestor_variables = std::vector<size_t>();
  another.unnormalised_log_weights = arma::colvec();
  another.normalised_log_weights = arma::colvec();
  another.previous_normalised_log_weights = arma::colvec();
  another.incremental_log_weights = arma::colvec();
  another.log_normalising_constant_ratio = 0.0;
  another.resampled_flag = false;
  another.schedule_parameters = Parameters();
  another.ess = 0.0;
}

/*
double& Particles::operator[](const size_t &i)
{
  return 1.0;
  //return(this->particles[i]);
}

double Particles::operator[](const size_t &i) const
{
  return 1.0;
  //return(this->particles[i]);
}
 */

MoveOutput* Particles::operator[](const size_t &i)
{
  return(this->particles[i]);
}

MoveOutput* Particles::operator[](const size_t &i) const
{
  return(this->particles[i]);
}

MoveOutput* Particles::back()
{
  return this->particles.back();
}

MoveOutput* Particles::back() const
{
  return this->particles.back();
}

size_t Particles::size() const
{
  return(this->particles.size());
}

/*
arma::colvec Particles::get_resampling_variables() const
{
  return this->resampling_variables;
}

void Particles::set_ancestor_variables(const std::vector<size_t> &ancestor_variables_in)
{
  this->ancestor_variables = ancestor_variables_in;
}
 */

void Particles::resample()
{
  this->ancestor_variables = stratified_resample(this->normalised_log_weights,
                                                 this->resampling_variables);
}

arma::rowvec Particles::get_output_lengths() const
{
  arma::rowvec output_lengths(this->size());
  
  for (size_t i=0; i<this->particles.size(); ++i)
  {
    output_lengths[i] = this->particles[i]->length();
  }
  
  return output_lengths;
}

void Particles::initialise_weights()
{
  this->log_normalising_constant_ratio = 0.0;
  this->previous_normalised_log_weights = arma::colvec(this->particles.size());
  this->previous_normalised_log_weights.fill(-log(double(this->particles.size())));
}

void Particles::update_weights(const arma::colvec &latest_unnormalised_log_incremental_weights)
{
  
  //arma::colvec latest_unnormalised_log_weights;
  if ( (latest_unnormalised_log_incremental_weights.size()>0) && (this->normalised_log_weights.size()>0) && (latest_unnormalised_log_incremental_weights.size()==this->normalised_log_weights.size()) )
    this->unnormalised_log_weights = this->previous_normalised_log_weights + latest_unnormalised_log_incremental_weights;
  else
    Rcpp::stop("SMCOutput::update_unnormalised_log_weights: weights have the wrong length.");
  
  //this->unnormalised_log_weights = this->previous_normalised_log_weights + latest_unnormalised_log_incremental_weights;
  
  this->incremental_log_weights = std::move(latest_unnormalised_log_incremental_weights);
  
  //this->unnormalised_log_weights = latest_unnormalised_log_weights;
  //this->log_normalising_constant_ratio = log_sum_exp(this->unnormalised_log_weights);
  
  // Set each particle to know what its probability was last time it was evaluated at the most recent target.
  /*
  for (std::vector< std::deque<Particle> >::iterator i = this->particles.begin();
       i != this->particles.end();
       ++i)
  {
    i->back().previous_target_evaluated = i->back().target_evaluated;
  }
  */
  
  //this->most_recent_log_ratio_estimator =
  
  //this->log_likelihood = this->log_likelihood + this->most_recent_log_ratio_estimator;
  
  //size_t num_to_pop_back = std::max<int>(0,unnormalised_log_weights.size()-lag);
  //for (size_t i=0; i<num_to_pop_back; ++i)
  //{
  //  this->unnormalised_log_weights.pop_back();
  //}
}

double Particles::calculate_log_normalising_constant()
{
  this->log_normalising_constant_ratio = log_sum_exp(this->unnormalised_log_weights);
  if (this->log_normalising_constant_ratio==-arma::datum::inf)
    this->ess = 0.0;
  else
    this->ess = exp(2.0*this->log_normalising_constant_ratio - log_sum_exp(2.0*this->unnormalised_log_weights));
  return this->log_normalising_constant_ratio;
}

void Particles::normalise_weights()
{
  this->normalised_log_weights = this->unnormalised_log_weights - this->log_normalising_constant_ratio;
}

arma::mat Particles::get_most_recent_matrix_particles(const std::vector<std::string> &variables) const
{
  if (this->particles.size()==0)
    return arma::mat(0,0);
  
  arma::rowvec first_vector = this->particles.front()->back().get_rowvec(variables);
  
  arma::mat output(this->particles.size(),first_vector.n_rows);
  output.row(0) = first_vector;
  for (size_t i=1; i<this->particles.size(); ++i)
  {
    output.row(i) = this->particles[i]->back().get_rowvec(variables);
  }
  return output;
}

void Particles::close_ofstreams()
{
  for (size_t i=0; i<this->particles.size(); ++i)
  {
    this->particles[i]->close_ofstreams();
  }
}

void Particles::set_previous_target_evaluated_to_target_evaluated()
{
  for (std::vector< MoveOutput* >::iterator i = this->particles.begin();
       i != this->particles.end();
       ++i)
  {
    (*i)->back().previous_target_evaluated = (*i)->back().target_evaluated;
    //(*i)->back().previous_target_gradients_of_log = (*i)->back().target_gradients_of_log;
  }
}

void Particles::subsample_set_previous_target_evaluated_to_target_evaluated()
{
  for (std::vector< MoveOutput* >::iterator i = this->particles.begin();
       i != this->particles.end();
       ++i)
  {
    (*i)->back().subsample_previous_target_evaluated = (*i)->back().subsample_target_evaluated;
    //(*i)->back().subsample_previous_target_gradients_of_log = (*i)->back().subsample_target_gradients_of_log;
  }
}

std::ostream& operator<<(std::ostream& os, const Particles &p)
{
  //std::vector< std::deque<Particle> >::const_iterator it;

  /*
  for (it=p.particles.begin();it!=p.particles.end();++it)
  {
    if (it==p.particles.begin())
      os << *it;
    else
      os << std::endl << *it;
  }
   */

  return os;
}
