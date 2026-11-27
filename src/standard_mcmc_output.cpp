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
  this->termination.reset();
  this->storage_decided = false;
  this->rich_storage = false;
  this->thin = 1;
  this->number_of_iterations_seen = 0;
}

StandardMCMCOutput::StandardMCMCOutput(MCMCTermination* termination_in)
{
  this->iteration_counter = 0;
  this->termination.reset(termination_in);
  this->storage_decided = false;
  this->rich_storage = false;
  this->thin = 1;
  this->number_of_iterations_seen = 0;
}

bool StandardMCMCOutput::reset_for_new_chain(const MCMC* mcmc_in, Particle &scratch)
{
  if ((this->get_mcmc() != mcmc_in) || this->records_adaptation || this->diagnostics)
    return false;
  this->iteration_counter = 0;
  if (this->storage_decided)
    this->store = SampleStore(); // (a fresh store: its layout and extra columns are decided again by the new chain)
  this->output.clear();
  this->storage_decided = false;
  this->rich_storage = false;
  this->number_of_iterations_seen = 0;
  this->adaptation_trace.clear();
  this->adaptation_final.clear();
  this->adaptation_columns = 0;
  this->algorithm_parameters.clear();
  if (this->termination != NULL)
    this->termination->set_parameters(this); // (its counter: this output's)
  scratch = std::move(this->current);
  return true;
}

StandardMCMCOutput::~StandardMCMCOutput()
{
  this->termination.reset();
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
  
  this->termination.reset();
  
  MoveOutput::operator=(another);
  this->make_copy(another);
}

void StandardMCMCOutput::make_copy(const StandardMCMCOutput &another)
{
  if (another.termination!=NULL)
  {
    this->termination.reset(another.termination->duplicate());
    // the copy's rule reads the copy's own iteration counter (it read the
    // original's: harmless while copies were never run again; a copied output
    // reused for a new chain then never stopped, C262)
    this->termination->set_parameters(this);
  }
  else
    this->termination.reset();
  
  this->iteration_counter = another.iteration_counter;
  this->algorithm_parameters = another.algorithm_parameters;
  this->output = another.output;
  this->current = another.current;
  this->store = another.store;
  this->storage_decided = another.storage_decided;
  this->rich_storage = another.rich_storage;
  this->thin = another.thin;
  this->number_of_iterations_seen = another.number_of_iterations_seen;
  this->records_adaptation = another.records_adaptation;
  this->adaptation_columns = another.adaptation_columns;
  this->adaptation_trace = another.adaptation_trace;
  this->adaptation_final = another.adaptation_final;
  this->diagnostics.reset(another.diagnostics ? new Diagnostics(*another.diagnostics) : NULL);
  this->signs_recorded = another.signs_recorded;
  this->signs = another.signs;
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

  // a stored iteration: with it, the adaptation state after this iteration
  if (this->records_adaptation)
  {
    size_t before = this->adaptation_trace.size();
    this->get_mcmc()->append_adaptation_state(NULL, &this->adaptation_trace, false, "");
    this->adaptation_columns = this->adaptation_trace.size() - before;
  }
  if (this->diagnostics)
    this->diagnostics->rows.insert(this->diagnostics->rows.end(), this->diagnostics->pending.begin(), this->diagnostics->pending.end());

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
    // (dev/plan.md 3.22 U9 R2) signed pseudo-marginal: each stored point's sign
    this->signs_recorded = (particle_in.factor_variables!=NULL) && particle_in.factor_variables->may_be_negative();
    this->storage_decided = true;
  }
  if (this->signs_recorded)
    this->signs.push_back(particle_in.factor_variables->sign());

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

bool StandardMCMCOutput::repoint_to(const Factors* factors)
{
  if ((this->current.factor_variables!=NULL) && !this->current.factor_variables->repoint_to(factors))
    return false;
  for (auto i = this->output.begin(); i != this->output.end(); ++i)
  {
    if ((i->factor_variables!=NULL) && !i->factor_variables->repoint_to(factors))
      return false;
  }
  return true;
}

void StandardMCMCOutput::start_recording_adaptation()
{
  this->records_adaptation = true;
  this->adaptation_trace.clear();
  this->adaptation_final.clear();
  this->adaptation_columns = 0;
}

void StandardMCMCOutput::start_recording_diagnostics(const std::vector<std::string> &names)
{
  this->diagnostics.reset(new Diagnostics());
  this->diagnostics->names = names;
  this->diagnostics->pending.assign(names.size(), 0.0);
  for (size_t c=0; c<names.size(); ++c)
  {
    const std::string &name = names[c];
    if ((name == "divergent") || ((name.size() > 10) && (name.compare(name.size() - 10, 10, ".divergent") == 0)))
      this->diagnostics->divergent_columns.push_back(c);
  }
}

void StandardMCMCOutput::note_diagnostics(const double* row, bool after_warm_up)
{
  Diagnostics &d = *this->diagnostics;
  std::copy(row, row + d.pending.size(), d.pending.begin());
  if (after_warm_up)
  {
    d.iterations_after_warm_up = d.iterations_after_warm_up + 1;
    // (an iteration of a scan is divergent if one of its moves was)
    bool divergent = false;
    for (size_t c : d.divergent_columns)
      divergent = divergent || ((row[c] != 0.0) && !std::isnan(row[c]));
    if (divergent)
      d.divergent_after_warm_up = d.divergent_after_warm_up + 1;
  }
}

size_t StandardMCMCOutput::fill_diagnostics_transposed(arma::mat &destination,
                                                       size_t first_column) const
{
  size_t columns = this->diagnostics_columns();
  if (columns == 0)
    return 0;
  if (destination.n_rows != columns)
    Rcpp::stop("StandardMCMCOutput: chains with different diagnostics cannot be written together.");
  size_t rows = this->diagnostics->rows.size()/columns;
  for (size_t r=0; r<rows; ++r)
    for (size_t c=0; c<columns; ++c)
      destination(c, first_column + r) = this->diagnostics->rows[r*columns + c];
  return rows;
}

size_t StandardMCMCOutput::fill_signs(std::vector<double> &destination) const
{
  destination.insert(destination.end(), this->signs.begin(), this->signs.end());
  return this->signs.size();
}

void StandardMCMCOutput::diagnostics_names(std::vector<std::string> &names) const
{
  names.clear();
  if (this->diagnostics)
    names = this->diagnostics->names;
}

void StandardMCMCOutput::count_divergences(size_t &divergent, size_t &iterations) const
{
  if (!this->diagnostics)
    return;
  divergent = divergent + this->diagnostics->divergent_after_warm_up;
  iterations = iterations + this->diagnostics->iterations_after_warm_up;
}

void StandardMCMCOutput::finish_adaptation()
{
  this->adaptation_final.clear();
  this->get_mcmc()->append_adaptation_state(NULL, &this->adaptation_final, true, "");
}

size_t StandardMCMCOutput::adaptation_trace_columns() const
{
  return this->records_adaptation ? this->adaptation_columns : 0;
}

size_t StandardMCMCOutput::adaptation_final_columns() const
{
  return this->records_adaptation ? this->adaptation_final.size() : 0;
}

size_t StandardMCMCOutput::fill_adaptation_trace_transposed(arma::mat &destination,
                                                            size_t first_column) const
{
  if (this->adaptation_columns==0)
    return 0;
  size_t rows = this->adaptation_trace.size()/this->adaptation_columns;
  if (destination.n_rows!=this->adaptation_columns)
    Rcpp::stop("StandardMCMCOutput: chains with different adaptation outputs cannot be written together.");
  for (size_t r=0; r<rows; ++r)
    for (size_t c=0; c<this->adaptation_columns; ++c)
      destination(c, first_column + r) = this->adaptation_trace[r*this->adaptation_columns + c];
  return rows;
}

void StandardMCMCOutput::fill_adaptation_final(arma::mat &destination,
                                               size_t column) const
{
  if (destination.n_rows!=this->adaptation_final.size())
    Rcpp::stop("StandardMCMCOutput: chains with different adaptation outputs cannot be written together.");
  for (size_t c=0; c<this->adaptation_final.size(); ++c)
    destination(c, column) = this->adaptation_final[c];
}

void StandardMCMCOutput::adaptation_names(std::vector<std::string> &trace_names,
                                          std::vector<std::string> &final_names) const
{
  // const: the names do not change the move
  MCMC* mcmc = const_cast<StandardMCMCOutput*>(this)->get_mcmc();
  mcmc->append_adaptation_state(&trace_names, NULL, false, "");
  mcmc->append_adaptation_state(&final_names, NULL, true, "");
}

}
