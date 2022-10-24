#include "pdmp_mcmc_output.h"

PDMPMCMCOutput::PDMPMCMCOutput()
  :MoveOutput()
{
}

PDMPMCMCOutput::~PDMPMCMCOutput()
{
  
}

//Copy constructor for the PDMPMCMCOutput class.
PDMPMCMCOutput::PDMPMCMCOutput(const PDMPMCMCOutput &another)
  :MoveOutput(another)
{
  this->make_copy(another);
}

void PDMPMCMCOutput::operator=(const PDMPMCMCOutput &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  MoveOutput::operator=(another);
  this->make_copy(another);
}

MoveOutput* PDMPMCMCOutput::duplicate(void)const
{
  return( new PDMPMCMCOutput(*this));
}

void PDMPMCMCOutput::make_copy(const PDMPMCMCOutput &another)
{
  this->dummy = another.dummy;
}

Particle& PDMPMCMCOutput::back()
{
  return dummy;
}

Particle PDMPMCMCOutput::back() const
{
  return dummy;
}

std::vector<Parameters> PDMPMCMCOutput::get_vector_of_parameters() const
{
  Rcpp::stop("PDMPMCMCOutput::get_vector_of_parameters() - not written yet.");
}
