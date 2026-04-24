#ifndef GAUSSIANNOISEPROPOSALKERNEL_H
#define GAUSSIANNOISEPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "proposal_kernel.h"
#include "particle.h"
#include "distributions.h"
#include "ilike_header.h"

namespace ilike
{
  /**
   * @file gaussian_noise_proposal_kernel.h
   * @brief Defines the GaussianNoiseProposalKernel class.
   *
   * A gaussian noise proposal kernel. Proposes new parameter values during MCMC or SMC moves using a gaussian noise distribution centred on the current state.
   *
   * @namespace ilike
   * @class GaussianNoiseProposalKernel
   * @brief A gaussian noise proposal kernel derived from ProposalKernel.
   */


class GaussianNoiseProposalKernel : public ProposalKernel
{
  
public:
  
  /**
   * @brief Default constructor for GaussianNoiseProposalKernel.
   */
  GaussianNoiseProposalKernel();
  /**
   * @brief Destructor for GaussianNoiseProposalKernel.
   */
  virtual ~GaussianNoiseProposalKernel();
  
  /**
   * @brief Copy constructor for GaussianNoiseProposalKernel.
   *
   * @param another The GaussianNoiseProposalKernel instance to copy from.
   */
  GaussianNoiseProposalKernel(const GaussianNoiseProposalKernel &another);
  
  /**
   * @brief Assignment operator for GaussianNoiseProposalKernel.
   *
   * @param another The GaussianNoiseProposalKernel instance to copy from.
   */
  void operator=(const GaussianNoiseProposalKernel &another);
  //Kernel* duplicate() const;
  /**
   * @brief Creates a deep copy and returns it as a gaussian_noise_proposal_kernel pointer.
   *
   * @return The result.
   */
  virtual GaussianNoiseProposalKernel* gaussian_noise_proposal_kernel_duplicate() const=0;
  
  /**
   * @brief Performs the can be evaluated operation.
   *
   * @return The result.
   */
  bool can_be_evaluated() const;
  
  /**
   * @brief Sets the data.
   *
   * @param data_in The data.
   */
  void set_data(Data* data_in);
  
protected:
  
  /**
   * @brief Copies the state of another GaussianNoiseProposalKernel into this object.
   *
   * @param another The GaussianNoiseProposalKernel instance to copy from.
   */
  void make_copy(const GaussianNoiseProposalKernel &another);
  
};
}

#endif
