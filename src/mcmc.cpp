#include "mcmc.h"
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
  this->termination = NULL;
  this->cut_index = NULL;
  this->index_with_cut = NULL;
  this->check_factor_index = false;
}

MCMC::MCMC(size_t number_of_iterations_in)
:Kernel()
{
  this->termination = new IterationsMCMCTermination(number_of_iterations_in);
  this->cut_index = NULL;
  this->index_with_cut = NULL;
  this->check_factor_index = false;
  //this->termination->set_parameters(this);
}

MCMC::MCMC(MCMCTermination* termination_in)
:Kernel()
{
  this->termination = termination_in;
  this->cut_index = NULL;
  this->index_with_cut = NULL;
  this->check_factor_index = false;
  //this->termination->set_parameters(this);
}

MCMC::~MCMC()
{
  if (this->termination!=NULL)
    delete this->termination;
  if (this->cut_index!=NULL)
    delete this->cut_index;
  if (this->index_with_cut!=NULL)
    delete this->index_with_cut;
  
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
  
  if (this->termination!=NULL)
    delete this->termination;
  
  //if (this->adaptor!=NULL)
  //  delete this->adaptor;
  
  if (this->cut_index!=NULL)
    delete this->cut_index;
  if (this->index_with_cut!=NULL)
    delete this->index_with_cut;
  Kernel::operator=(another);
  this->make_copy(another);
}

void MCMC::make_copy(const MCMC &another)
{
  if (another.termination!=NULL)
    this->termination = another.termination->duplicate();
  else
    this->termination = NULL;
  
  //if (another.adaptor!=NULL)
  //  this->adaptor = another.adaptor->duplicate();
  if (another.cut_index!=NULL)
    this->cut_index = another.cut_index->duplicate();
  else
    this->cut_index = NULL;
  if (another.index_with_cut!=NULL)
    this->index_with_cut = another.index_with_cut->duplicate();
  else
    this->index_with_cut = NULL;
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
  if (this->cut_index!=NULL)
    delete this->cut_index;
  this->cut_index = (cut_index_in!=NULL) ? cut_index_in->duplicate() : NULL;
}

void MCMC::update_index_with_cut(const Index* move_index)
{
  if (this->index_with_cut!=NULL)
  {
    delete this->index_with_cut;
    this->index_with_cut = NULL;
  }
  if ( (this->cut_index==NULL) || (move_index==NULL) )
    return;
  std::vector<size_t> both(move_index->begin(), move_index->end());
  for (auto i=this->cut_index->begin(); i!=this->cut_index->end(); ++i)
  {
    if (std::find(both.begin(), both.end(), *i)==both.end())
      both.push_back(*i);
  }
  this->index_with_cut = new VectorIndex(both);
}

const Index* MCMC::get_new_outputs_index(const Index* move_index) const
{
  if (this->index_with_cut!=NULL)
    return this->index_with_cut;
  else
    return move_index;
}

void MCMC::refresh_cut_factors(Particle &accepted_particle) const
{
  if (this->cut_index!=NULL)
    accepted_particle.evaluate_likelihoods(this->cut_index);
}

MoveOutput* MCMC::run(RandomNumberGenerator &rng,
                      const Particle &particle) const

{
  StandardMCMCOutput* mcmc_output = this->initialise_mcmc_output();//new StandardMCMCOutput();
  
  Particle current_particle = particle;
  Particle* old_pointer = &current_particle;
  while (!mcmc_output->terminate())
  {
    current_particle = mcmc_output->move(rng,
                                         current_particle);
    mcmc_output->increment_counter();
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
