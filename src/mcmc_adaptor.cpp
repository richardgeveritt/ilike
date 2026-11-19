#include "mcmc_adaptor.h"
#include "proposal_kernel.h"

namespace ilike
{

MCMCAdaptor::MCMCAdaptor()
{
  this->proposal = NULL;
}

MCMCAdaptor::~MCMCAdaptor()
{
}

MCMCAdaptor::MCMCAdaptor(const MCMCAdaptor &another)
{
  this->make_copy(another);
}

void MCMCAdaptor::operator=(const MCMCAdaptor &another)
{
  if(this == &another)
    return;
  this->make_copy(another);
}

void MCMCAdaptor::make_copy(const MCMCAdaptor &another)
{
  // the copy belongs to a copy of the kernel, which points it at itself
  // (ProposalKernel::make_copy); until then it points at the original's
  this->proposal = another.proposal;
  this->until.reset((another.until != NULL) ? another.until->duplicate() : NULL);
  this->stopped = another.stopped;
  this->steps = another.steps;
  this->stopped_at_step = another.stopped_at_step;
}

void MCMCAdaptor::set_proposal(ProposalKernel* proposal_in)
{
  this->proposal = proposal_in;
}

void MCMCAdaptor::prepare_for_move()
{
}

void MCMCAdaptor::set_until(AdaptationCriterion* until_in)
{
  this->until.reset(until_in);
}

void MCMCAdaptor::adapt_unless_stopped(const AdaptationStep &step)
{
  if (this->stopped)
    return;
  this->adapt(step);
  this->steps = this->steps + 1;
  if ((this->until != NULL) && this->until->met(*this, step.acceptance_probability, this->steps))
  {
    this->stopped = true;
    this->stopped_at_step = this->steps;
  }
}

bool MCMCAdaptor::watched(const std::string &name, std::vector<double> &values) const
{
  return false;
}

std::vector<std::string> MCMCAdaptor::watchable() const
{
  return std::vector<std::string>();
}

std::vector<std::string> MCMCAdaptor::default_watch() const
{
  return this->watchable();
}

bool MCMCAdaptor::acceptance_target(double &target) const
{
  return false;
}

std::string MCMCAdaptor::state_label() const
{
  return "the proposal";
}

void MCMCAdaptor::mcmc_adapt(const Particle &latest_particle,
                             size_t iteration_counter)
{
  Rcpp::stop("MCMC adaptation is available for Metropolis (m_proposal) and Metropolis-Hastings (mh_proposal) moves only.");
}

}
