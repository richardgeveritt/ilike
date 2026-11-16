#include "mcmc.h"
#include "chunk_parallel.h"
#include "mcmc_termination.h"
#include "iterations_mcmc_termination.h"
#include "move_output.h"
#include "standard_mcmc_output.h"
#include "mcmc_adaptor.h"
#include "stochastic_scan_mcmc.h"
#include "deterministic_scan_mcmc.h"
#include "index.h"
#include "vector_index.h"
#include "factor_variables.h"
#include <algorithm>

namespace ilike
{
MCMC::MCMC()
:Kernel()
{
  this->termination.reset();
  this->cut_index.reset();
  this->index_with_cut.reset();
  this->check_factor_index = false;
}

MCMC::MCMC(size_t number_of_iterations_in)
:Kernel()
{
  this->termination.reset(new IterationsMCMCTermination(number_of_iterations_in));
  this->cut_index.reset();
  this->index_with_cut.reset();
  this->check_factor_index = false;
  //this->termination->set_parameters(this);
}

MCMC::MCMC(MCMCTermination* termination_in)
:Kernel()
{
  this->termination.reset(termination_in);
  this->cut_index.reset();
  this->index_with_cut.reset();
  this->check_factor_index = false;
  //this->termination->set_parameters(this);
}

MCMC::~MCMC()
{
  this->termination.reset();
  this->cut_index.reset();
  this->index_with_cut.reset();
  
  //if (this->adaptor!=NULL)
  //  delete this->adaptor;
}

MCMC::MCMC(const MCMC &another)
:Kernel(another)
{
  this->make_copy(another);
}

void MCMC::operator=(const MCMC &another)
{
  if(this == &another)
    return;
  
  this->termination.reset();
  
  //if (this->adaptor!=NULL)
  //  delete this->adaptor;
  
  this->cut_index.reset();
  this->index_with_cut.reset();
  Kernel::operator=(another);
  this->make_copy(another);
}

void MCMC::make_copy(const MCMC &another)
{
  this->thin = another.thin;
  if (another.termination!=NULL)
    this->termination.reset(another.termination->duplicate());
  else
    this->termination.reset();
  
  //if (another.adaptor!=NULL)
  //  this->adaptor = another.adaptor->duplicate();
  if (another.cut_index!=NULL)
    this->cut_index.reset(another.cut_index->duplicate());
  else
    this->cut_index.reset();
  if (another.index_with_cut!=NULL)
    this->index_with_cut.reset(another.index_with_cut->duplicate());
  else
    this->index_with_cut.reset();
  this->check_factor_index = another.check_factor_index;
  this->move_label = another.move_label;
}

bool MCMC::reads_cached_factors() const
{
  return true;
}

void MCMC::set_check_factor_index(bool check, const std::string &label)
{
  this->check_factor_index = check;
  this->move_label = label;
}

void MCMC::check_carried_over_factors(const Particle &accepted_particle,
                                      const Index* move_index) const
{
  if ( (!this->check_factor_index) || (accepted_particle.factor_variables==NULL) )
    return;
  std::vector<size_t> changed = accepted_particle.factor_variables->find_changed_exact_factors(this->get_new_outputs_index(move_index));
  if (changed.size()>0)
  {
    std::string factors;
    for (size_t i=0; i<changed.size(); ++i)
      factors = factors + (i>0 ? ", " : "") + std::to_string(changed[i]+1);
    Rcpp::stop("check_factor_index: in " + this->move_label + ", factor(s) " + factors +
               " changed value when the move changed the parameters, so they do not cancel in its acceptance ratio. Add them to the factor index, or to the cut factor index to deliberately leave them out.");
  }
}

void MCMC::set_cut_index(const Index* cut_index_in)
{
  Rcpp::stop("MCMC::set_cut_index - cut factors are only available for Metropolis, Metropolis-Hastings and unadjusted moves.");
}

void MCMC::store_cut_index(const Index* cut_index_in)
{
  this->cut_index.reset();
  this->cut_index.reset((cut_index_in!=NULL) ? cut_index_in->duplicate() : NULL);
}

void MCMC::update_index_with_cut(const Index* move_index)
{
  if (this->index_with_cut!=NULL)
  {
    this->index_with_cut.reset();
    this->index_with_cut.reset();
  }
  if ( (this->cut_index==NULL) || (move_index==NULL) )
    return;
  std::vector<size_t> both(move_index->begin(), move_index->end());
  for (auto i=this->cut_index->begin(); i!=this->cut_index->end(); ++i)
  {
    if (std::find(both.begin(), both.end(), *i)==both.end())
      both.push_back(*i);
  }
  this->index_with_cut.reset(new VectorIndex(both));
}

const Index* MCMC::get_new_outputs_index(const Index* move_index) const
{
  if (this->index_with_cut!=NULL)
    return this->index_with_cut.get();
  else
    return move_index;
}

void MCMC::refresh_cut_factors(Particle &accepted_particle) const
{
  if (this->cut_index!=NULL)
    accepted_particle.evaluate_likelihoods(this->cut_index.get());
}

MoveOutput* MCMC::run(RandomNumberGenerator &rng,
                      const Particle &particle) const

{
  StandardMCMCOutput* mcmc_output = this->initialise_mcmc_output();//new StandardMCMCOutput();
  
  mcmc_output->set_thin(this->thin);
  bool adapting = this->adapts();
  if (adapting)
    mcmc_output->start_recording_adaptation();
  Particle current_particle = particle;
  // proposals are built here, reusing its objects (Stage 3 (a))
  Particle scratch;
  while (!mcmc_output->terminate())
  {
    // In place: a rejected proposal no longer copies the current state.
    if (adapting)
      mcmc_output->get_mcmc()->adaptive_move_in_place(rng, current_particle, scratch);
    else
      mcmc_output->get_mcmc()->move_in_place(rng, current_particle, scratch);
    mcmc_output->increment_counter();
    if ((*mcmc_output->get_iteration_counter_pointer() % 1024) == 0)
      check_interrupt(); // in long chains (dev/plan.md 3.11)
    current_particle.tell_factors_to_forget_they_were_already_written_to_file();
    // Store the iteration compactly (its values and factor records), not a
    // copy of the particle; the state itself is handed over at the end.
    mcmc_output->record(current_particle);
  }
  if (adapting)
    mcmc_output->finish_adaptation();
  mcmc_output->set_current(std::move(current_particle));
  mcmc_output->reset_counter();
  return mcmc_output;
}

size_t MCMC::known_number_of_iterations() const
{
  if (this->termination==NULL)
    return 0;
  return this->termination->get_known_number_of_iterations();
}

void MCMC::set_thin(size_t thin_in)
{
  if (thin_in < 1)
    Rcpp::stop("MCMC: thin must be at least 1.");
  this->thin = thin_in;
}

void MCMC::move_in_place(RandomNumberGenerator &rng,
                         Particle &particle,
                         Particle &scratch) const
{
  this->move_in_place(rng, particle);
}

void MCMC::append_adaptation_state(std::vector<std::string> *names,
                                   std::vector<double> *values,
                                   bool final,
                                   const std::string &prefix) const
{
}

void MCMC::adaptive_move_in_place(RandomNumberGenerator &rng,
                                  Particle &particle,
                                  Particle &scratch)
{
  this->move_in_place(rng, particle, scratch);
  this->mcmc_adapt(particle, 0);
}

void MCMC::move_in_place(RandomNumberGenerator &rng,
                         Particle &particle) const
{
  particle = this->move(rng, particle);
}

/*
 MoveOutput* MCMC::run(RandomNumberGenerator &rng,
 Particle &particle,
 const Parameters &conditioned_on_parameters)
 
 {
 StandardMCMCOutput* mcmc_output = new StandardMCMCOutput();
 Particle current_particle = particle;
 while (!this->termination->terminate())
 {
 current_particle = this->move(rng,
 current_particle,
 conditioned_on_parameters);
 this->iteration_counter = this->iteration_counter + 1;
 this->mcmc_adapt(current_particle,
 this->iteration_counter);
 mcmc_output->push_back(current_particle);
 }
 return mcmc_output;
 }
 */

MoveOutput* MCMC::subsample_run(RandomNumberGenerator &rng,
                                const Particle &particle) const

{
  StandardMCMCOutput* mcmc_output = this->initialise_mcmc_output();
  
  Particle current_particle = particle;
  Particle* old_pointer = &current_particle;
  while (!mcmc_output->terminate())
  {
    current_particle = this->subsample_move(rng,
                                            current_particle);
    mcmc_output->increment_counter();
    if ((*mcmc_output->get_iteration_counter_pointer() % 1024) == 0)
      check_interrupt(); // in long chains (dev/plan.md 3.11)
    if (&current_particle==old_pointer)
    {
      current_particle.tell_factors_to_forget_they_were_already_written_to_file();
    }
    mcmc_output->push_back(current_particle);
    mcmc_output->mcmc_adapt();
  }
  mcmc_output->reset_counter();
  return mcmc_output;
}

void MCMC::mcmc_adapt(const Particle &current_particle,
                      size_t iteration_counter)
{
  this->specific_mcmc_adapt(current_particle,
                            iteration_counter);
  
  //Rcpp::stop("MCMC::mcmc_adapt - need to clear adaptation info.");
  //current_particle.erase_mcmc_adaptation_info();
}

MCMCTermination* MCMC::get_duplicated_termination() const
{
  return this->termination->duplicate();
}

/*
 MoveOutput* MCMC::subsample_run(RandomNumberGenerator &rng,
 Particle &particle,
 const Parameters &conditioned_on_parameters)
 
 {
 StandardMCMCOutput* mcmc_output = new StandardMCMCOutput();
 Particle current_particle = particle;
 while (!this->termination->terminate())
 {
 current_particle = this->subsample_move(rng,
 current_particle,
 conditioned_on_parameters);
 this->mcmc_adapt(current_particle,
 this->iteration_counter);
 this->iteration_counter = this->iteration_counter + 1;
 mcmc_output->push_back(current_particle);
 }
 return mcmc_output;
 }
 */

/*
 EnsembleMember MCMC::run(RandomNumberGenerator &rng,
 EnsembleMember &particle)
 
 {
 StandardMCMCOutput* mcmc_output = new StandardMCMCOutput();
 EnsembleMember current_particle = particle;
 while (!this->termination->terminate())
 {
 current_particle = this->move(rng,
 current_particle);
 this->iteration_counter = this->iteration_counter + 1;
 this->mcmc_adapt(&current_particle,
 this->iteration_counter);
 }
 return current_particle;
 }
 
 EnsembleMember MCMC::run(RandomNumberGenerator &rng,
 EnsembleMember &particle,
 const Parameters &conditioned_on_parameters)
 
 {
 StandardMCMCOutput* mcmc_output = new StandardMCMCOutput();
 EnsembleMember current_particle = particle;
 while (!this->termination->terminate())
 {
 current_particle = this->move(rng,
 current_particle,
 conditioned_on_parameters);
 this->iteration_counter = this->iteration_counter + 1;
 this->mcmc_adapt(&current_particle,
 this->iteration_counter);
 }
 return current_particle;
 }
 
 EnsembleMember MCMC::subsample_run(RandomNumberGenerator &rng,
 EnsembleMember &particle,
 const Parameters &conditioned_on_parameters)
 
 {
 StandardMCMCOutput* mcmc_output = new StandardMCMCOutput();
 EnsembleMember current_particle = particle;
 while (!this->termination->terminate())
 {
 current_particle = this->subsample_move(rng,
 current_particle,
 conditioned_on_parameters);
 this->mcmc_adapt(&current_particle,
 this->iteration_counter);
 this->iteration_counter = this->iteration_counter + 1;
 }
 return current_particle;
 }
 */
}
