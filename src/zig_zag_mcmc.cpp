#include <cmath>
#include "distributions.h"
#include "utils.h"
#include "zig_zag_mcmc.h"
#include "pdmp_mcmc_output.h"
#include "metropolis_standard_mcmc_output.h" // will be changed

namespace ilike
{
ZigZagMCMC::ZigZagMCMC()
:MCMC()
{
}

ZigZagMCMC::ZigZagMCMC(size_t number_of_iterations_in,
                       ProposalKernel* proposal_in)
:MCMC(number_of_iterations_in)
{
  this->proposal = proposal_in;
}

ZigZagMCMC::ZigZagMCMC(size_t number_of_iterations_in,
                       const std::vector<Parameters> &initial_points_in,
                       const Parameters &proposal_variances_in)
:MCMC(number_of_iterations_in)
{
  // default to Gaussian random walk
  //this->proposal = ProposalKernel(EvaluateLogMCMCProposalPtr proposal_evaluate_in,
  //                                SimulateMCMCProposalPtr proposal_simulate_in,
  //                                proposal_variances_in);
}

ZigZagMCMC::~ZigZagMCMC()
{
  if (this->proposal!=NULL)
    delete this->proposal;
}

//Copy constructor for the ZigZagMCMC class.
ZigZagMCMC::ZigZagMCMC(const ZigZagMCMC &another)
:MCMC(another)
{
  this->make_copy(another);
}

void ZigZagMCMC::operator=(const ZigZagMCMC &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  if (this->proposal!=NULL)
    delete this->proposal;
  
  MCMC::operator=(another);
  this->make_copy(another);
}

Kernel* ZigZagMCMC::duplicate() const
{
  return( new ZigZagMCMC(*this));
}

MCMC* ZigZagMCMC::mcmc_duplicate() const
{
  return( new ZigZagMCMC(*this));
}

void ZigZagMCMC::make_copy(const ZigZagMCMC &another)
{
  this->proposal = another.proposal;
  if (another.proposal!=NULL)
    this->proposal = another.proposal->proposal_kernel_duplicate();
  else
    this->proposal = NULL;
}

Particle ZigZagMCMC::move(RandomNumberGenerator &rng,
                          const Particle &particle) const
{
  Rcpp::stop("ZigZagMCMC::move - not yet implemented.");
}

/*
 Particle ZigZagMCMC::move(RandomNumberGenerator &rng,
 Particle &particle,
 const Parameters &conditioned_on_parameters) const
 {
 Rcpp::stop("ZigZagMCMC::move - not yet implemented.");
 }
 */

Particle ZigZagMCMC::subsample_move(RandomNumberGenerator &rng,
                                    const Particle &particle) const
{
  Rcpp::stop("ZigZagMCMC::move - not yet implemented.");
}

/*
 Particle ZigZagMCMC::subsample_move(RandomNumberGenerator &rng,
 Particle &particle,
 const Parameters &conditioned_on_parameters) const
 {
 Rcpp::stop("ZigZagMCMC::move - not yet implemented.");
 }
 */

/*
 EnsembleMember ZigZagMCMC::move(RandomNumberGenerator &rng,
 const Index* index,
 EnsembleMember &particle) const
 {
 Rcpp::stop("ZigZagMCMC::move - not yet implemented.");
 }
 
 EnsembleMember ZigZagMCMC::move(RandomNumberGenerator &rng,
 const Index* index,
 EnsembleMember &particle,
 const Parameters &conditioned_on_parameters) const
 {
 Rcpp::stop("ZigZagMCMC::move - not yet implemented.");
 }
 
 EnsembleMember ZigZagMCMC::subsample_move(RandomNumberGenerator &rng,
 const Index* index,
 EnsembleMember &particle,
 const Parameters &conditioned_on_parameters) const
 {
 Rcpp::stop("ZigZagMCMC::move - not yet implemented.");
 }
 */

MoveOutput* ZigZagMCMC::run(RandomNumberGenerator &rng,
                            const Particle &particle) const
{
  // run PDMP
  return new PDMPMCMCOutput();
}

/*
 MoveOutput* ZigZagMCMC::run(RandomNumberGenerator &rng,
 Particle &particle,
 const Parameters &conditioned_on_parameters) const
 {
 // run PDMP
 return new PDMPMCMCOutput();
 }
 */

MoveOutput* ZigZagMCMC::subsample_run(RandomNumberGenerator &rng,
                                      const Particle &particle) const
{
  // run PDMP
  return new PDMPMCMCOutput();
}

/*
 MoveOutput* ZigZagMCMC::subsample_run(RandomNumberGenerator &rng,
 Particle &particle,
 const Parameters &conditioned_on_parameters) const
 {
 // run PDMP
 return new PDMPMCMCOutput();
 }
 */

void ZigZagMCMC::smc_adapt(SMCOutput* current_state)
{
  proposal->smc_adapt(current_state);
}

void ZigZagMCMC::ensemble_adapt(EnsembleKalmanOutput* current_state)
{
  proposal->ensemble_adapt(current_state);
}

void ZigZagMCMC::specific_mcmc_adapt(const Particle &current_particle,
                                     size_t iteration_counter)
{
  proposal->mcmc_adapt(current_particle,
                       iteration_counter);
}

void ZigZagMCMC::set_index(Index* index_in)
{
  
}

void ZigZagMCMC::set_index_if_null(Index* index_in)
{
  
}

void ZigZagMCMC::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  this->proposal->set_proposal_parameters(proposal_parameters_in);
}

std::vector<const ProposalKernel*> ZigZagMCMC::get_proposals() const
{
  return std::vector<const ProposalKernel*>();
}

StandardMCMCOutput* ZigZagMCMC::initialise_mcmc_output() const
{
  return new MetropolisStandardMCMCOutput();
}
}
