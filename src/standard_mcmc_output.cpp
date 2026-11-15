#include "standard_mcmc_output.h"
#include "particle.h"
#include "factor_variables.h"
#include "ensemble_factor_variables.h"
#include "transform.h"
#include "mcmc_termination.h"
#include "mcmc.h"

namespace ilike
{
StandardMCMCOutput::StandardMCMCOutput()
:MoveOutput()
{
  this->iteration_counter = 0;
  this->termination = NULL;
  this->storage_decided = false;
  this->rich_storage = false;
  this->thin = 1;
  this->number_of_iterations_seen = 0;
}

StandardMCMCOutput::StandardMCMCOutput(MCMCTermination* termination_in)
{
  this->iteration_counter = 0;
  this->termination = termination_in;
  this->storage_decided = false;
  this->rich_storage = false;
  this->thin = 1;
  this->number_of_iterations_seen = 0;
}

//StandardMCMCOutput::StandardMCMCOutput(const Parameters &parameters_in)
//{
//  this->output = std::deque<Particle>(1);
//  this->output[0] = Particle(parameters_in);
//}

StandardMCMCOutput::~StandardMCMCOutput()
{
  if (this->termination!=NULL)
    delete this->termination;
}

//Copy constructor for the StandardMCMCOutput class.
StandardMCMCOutput::StandardMCMCOutput(const StandardMCMCOutput &another)
:MoveOutput(another)
{
  this->make_copy(another);
}

void StandardMCMCOutput::operator=(const StandardMCMCOutput &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  if (this->termination!=NULL)
    delete this->termination;
  
  MoveOutput::operator=(another);
  this->make_copy(another);
}

void StandardMCMCOutput::make_copy(const StandardMCMCOutput &another)
{
  if (another.termination!=NULL)
    this->termination = another.termination->duplicate();
  else
    this->termination = NULL;
  
  this->iteration_counter = another.iteration_counter;
  this->algorithm_parameters = another.algorithm_parameters;
  this->output = another.output;
  this->current = another.current;
  this->store = another.store;
  this->storage_decided = another.storage_decided;
  this->rich_storage = another.rich_storage;
  this->thin = another.thin;
  this->number_of_iterations_seen = another.number_of_iterations_seen;
}

void StandardMCMCOutput::push_back(const Particle &particle_in)
{
  this->record(particle_in);
  this->current = particle_in;
}

void StandardMCMCOutput::record(const Particle &particle_in)
{
  this->number_of_iterations_seen = this->number_of_iterations_seen + 1;
  if ((this->thin > 1) && ((this->number_of_iterations_seen % this->thin) != 0))
  {
      return;
  }

  if (!this->storage_decided)
  {
    // Compact storage needs every factor output to have a compact record,
    // and no ensemble factors (see the members in the header).
    this->rich_storage = (particle_in.ensemble_factor_variables!=NULL) ||
      ((particle_in.factor_variables!=NULL) && !particle_in.factor_variables->has_compact_values());
    if (!this->rich_storage)
    {
      size_t number_of_values = 0;
      if (particle_in.factor_variables!=NULL)
        number_of_values = particle_in.factor_variables->number_of_compact_values();
      this->store.set_number_of_extra_columns(number_of_values);
      this->compact_scratch.assign(number_of_values, 0.0);
      size_t known = (this->termination!=NULL) ? this->termination->get_known_number_of_iterations() : 0;
      if (known > 0)
        this->store.reserve_after_first(known/this->thin + 1);
    }
    this->storage_decided = true;
  }

  if (this->rich_storage)
  {
    this->output.push_back(particle_in);
  }
  else
  {
    if (this->compact_scratch.size() > 0)
    {
      particle_in.factor_variables->compact_values(this->compact_scratch.data());
      this->store.append(particle_in.parameters, this->compact_scratch.data());
    }
    else
    {
      this->store.append(particle_in.parameters);
    }
  }
}

void StandardMCMCOutput::set_current(Particle &&particle_in)
{
  this->current = std::move(particle_in);
}

void StandardMCMCOutput::set_thin(size_t thin_in)
{
  if (thin_in < 1)
    Rcpp::stop("MCMC: thin must be at least 1.");
  this->thin = thin_in;
}

Particle& StandardMCMCOutput::back()
{
  return this->current;
}

Particle StandardMCMCOutput::back() const
{
  return this->current;
}

std::vector<Parameters> StandardMCMCOutput::get_vector_of_parameters() const
{
  std::vector<Parameters> local_output;
  if (this->rich_storage)
  {
    local_output.reserve(this->output.size());
    for (auto i=this->output.begin(); i!=this->output.end(); ++i)
      local_output.push_back(i->parameters);
  }
  else
  {
    size_t n = this->store.size();
    local_output.reserve(n);
    for (size_t i=0; i<n; ++i)
      local_output.push_back(this->store.parameters_at(i));
  }
  return local_output;
}

void StandardMCMCOutput::write_vector_points(const std::vector<std::string> &variables,
                                             std::ofstream &file_stream,
                                             std::shared_ptr<Transform> transform) const
{
  if(file_stream.is_open())
  {
    std::vector<Parameters> points = this->get_vector_of_parameters();
    file_stream.precision(std::numeric_limits<double>::max_digits10);
    for (auto i = points.begin(); i!=points.end(); ++i)
    {
      if (transform==NULL)
        i->get_rowvec(variables).raw_print(file_stream);
      else
        transform->inverse_transform(*i).get_rowvec(variables).raw_print(file_stream);
    }
  }
}

void StandardMCMCOutput::write_any_points(const std::vector<std::string> &variables,
                                          std::ofstream &file_stream) const
{
  //if(file_stream.is_open())
  //{
  //  for (auto i = this->output.begin();
  //       i!=this->output.end();
  //       ++i)
  //  {
  //    file_stream << i->parameters. << std::endl;
  //  }
  //}
}

size_t StandardMCMCOutput::fill_vector_points_transposed(const std::vector<std::string> &variables,
                                                         arma::mat &destination,
                                                         size_t first_column) const
{
  if (this->rich_storage)
    return MoveOutput::fill_vector_points_transposed(variables, destination, first_column);
  size_t n = this->store.size();
  if (n > 0)
    this->store.fill_transposed(variables, destination, first_column);
  return n;
}

arma::mat StandardMCMCOutput::get_matrix_of_vector_points(const std::vector<std::string> &variables,
                                                           std::shared_ptr<Transform> transform) const
{
  if (this->rich_storage)
  {
    if (this->output.empty())
      return arma::mat();
    arma::rowvec first_row;
    if (transform == NULL)
      first_row = this->output.front().get_rowvec(variables);
    else
      first_row = transform->inverse_transform(this->output.front().parameters).get_rowvec(variables);
    arma::mat result(this->output.size(), first_row.n_elem);
    size_t row = 0;
    for (auto i = this->output.begin(); i != this->output.end(); ++i, ++row)
    {
      if (transform == NULL)
        result.row(row) = i->get_rowvec(variables);
      else
        result.row(row) = transform->inverse_transform(i->parameters).get_rowvec(variables);
    }
    return result;
  }

  size_t n = this->store.size();
  if (n == 0)
    return arma::mat();
  if (transform == NULL)
  {
    return this->store.matrix_of(variables);
  }
  std::vector<Parameters> points = this->get_vector_of_parameters();
  arma::rowvec first_row = transform->inverse_transform(points[0]).get_rowvec(variables);
  arma::mat result(n, first_row.n_elem);
  for (size_t i=0; i<n; ++i)
    result.row(i) = transform->inverse_transform(points[i]).get_rowvec(variables);
  return result;
}

void StandardMCMCOutput::write_factors(const std::string &directory_name,
                                       const std::string &index) const
{
  if (this->rich_storage)
  {
    for (auto i = this->output.begin();
         i!=this->output.end();
         ++i)
    {
      if (i->factor_variables!=NULL)
      {
        i->factor_variables->write_to_file(directory_name,
                                           index);
      }
    }
    return;
  }

  if ((this->current.factor_variables==NULL) || (this->store.get_number_of_extra_columns()==0))
    return;
  size_t n = this->store.size();
  // Every stored iteration is written with the values it had when it was
  // stored (the current state's factor outputs are only used to reach the
  // estimators that write them).
  for (size_t i=0; i<n; ++i)
    this->current.factor_variables->write_compact_values(directory_name, index, this->store.extras_of(i));
}

void StandardMCMCOutput::write_ensemble_factors(const std::string &directory_name,
                                                const std::string &index) const
{
  // Compact storage is only used without ensemble factor variables.
  for (auto i = this->output.begin();
       i!=this->output.end();
       ++i)
  {
    if (i->ensemble_factor_variables!=NULL)
    {
      i->ensemble_factor_variables->write_to_file(directory_name,
                                                  index);
    }
  }
}

size_t StandardMCMCOutput::length() const
{
  if (this->rich_storage)
    return this->output.size();
  return this->store.size();
}

void StandardMCMCOutput::close_ofstreams()
{
  if (this->rich_storage)
  {
    for (auto i = this->output.begin();
         i!=this->output.end();
         ++i)
    {
      if (i->factor_variables!=NULL)
      {
        i->factor_variables->close_ofstreams();
      }
      if (i->ensemble_factor_variables!=NULL)
      {
        i->ensemble_factor_variables->close_ofstreams();
      }
    }
    return;
  }
  // the stored records were written through the current state's factor
  // outputs (same estimators), so closing those flushes them
  if (this->current.factor_variables!=NULL)
    this->current.factor_variables->close_ofstreams();
}

size_t* StandardMCMCOutput::get_iteration_counter_pointer()
{
  return &this->iteration_counter;
}

void StandardMCMCOutput::increment_counter()
{
  this->iteration_counter = this->iteration_counter + 1;
}

void StandardMCMCOutput::reset_counter()
{
  this->iteration_counter = 0;
}

void StandardMCMCOutput::mcmc_adapt()
{
  this->get_mcmc()->mcmc_adapt(this->current,
                               this->iteration_counter);
}

void StandardMCMCOutput::mcmc_adapt(const Particle &current_particle)
{
  this->get_mcmc()->mcmc_adapt(current_particle,
                               this->iteration_counter);
}

bool StandardMCMCOutput::terminate() const
{
  return this->termination->terminate();
}

Particle StandardMCMCOutput::move(RandomNumberGenerator &rng,
                                  const Particle &particle) const
{
  return this->get_mcmc()->move(rng, particle);
}

Particle StandardMCMCOutput::subsample_move(RandomNumberGenerator &rng,
                                            const Particle &particle) const
{
  return this->get_mcmc()->subsample_move(rng, particle);
}

Parameters StandardMCMCOutput::get_current_algorithm_parameters() const
{
  return this->algorithm_parameters.back();
}
}
