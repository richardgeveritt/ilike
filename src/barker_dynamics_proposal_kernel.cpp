#include <math.h>
#include <iterator>
#include "barker_dynamics_proposal_kernel.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"
#include "distributions.h"
#include "gradient_estimator.h"
#include "independent_proposal_kernel.h"
#include "gradient_estimator_output.h"
#include "gaussian_independent_proposal_kernel.h"

namespace ilike
{
  BarkerDynamicsProposalKernel::BarkerDynamicsProposalKernel()
      : ProposalKernel()
  {
    this->gradient_estimator.reset();
    this->proposal_simulate.reset();
    this->index.reset();
  }

  BarkerDynamicsProposalKernel::~BarkerDynamicsProposalKernel()
  {
    this->gradient_estimator.reset();

    this->proposal_simulate.reset();

    this->index.reset();
  }

  BarkerDynamicsProposalKernel::BarkerDynamicsProposalKernel(const std::string &variable_name_in,
                                                             const arma::mat &covariance_in,
                                                             GradientEstimator *gradient_estimator_in)
      : ProposalKernel()
  {
    this->gradient_estimator.reset(gradient_estimator_in);
    this->gradient_estimator->set_proposal(this);
    this->index.reset();
    arma::colvec mean = arma::zeros<arma::colvec>(covariance_in.n_rows);
    // Create an identity matrix for the covariance
    arma::mat identity = arma::eye<arma::mat>(covariance_in.n_rows, covariance_in.n_rows);
    this->proposal_simulate.reset(new GaussianIndependentProposalKernel(variable_name_in,
                                                                    mean,
                                                                    identity));

    this->proposal_info[variable_name_in] = GaussianProposalInfo(covariance_in);
  }

  BarkerDynamicsProposalKernel::BarkerDynamicsProposalKernel(const std::string &variable_name_in,
                                                             const arma::mat &covariance_in,
                                                             double scale_in,
                                                             GradientEstimator *gradient_estimator_in)
      : ProposalKernel()
  {
    this->gradient_estimator.reset(gradient_estimator_in);
    this->gradient_estimator->set_proposal(this);
    this->index.reset();
    arma::colvec mean = arma::zeros<arma::colvec>(covariance_in.n_rows);
    // Create an identity matrix for the covariance
    arma::mat identity = arma::eye<arma::mat>(covariance_in.n_rows, covariance_in.n_rows);
    this->proposal_simulate.reset(new GaussianIndependentProposalKernel(variable_name_in,
                                                                    mean,
                                                                    identity));

    this->proposal_info[variable_name_in] = GaussianProposalInfo(covariance_in,
                                                                 scale_in);
  }

  BarkerDynamicsProposalKernel::BarkerDynamicsProposalKernel(const std::vector<std::string> &variable_names_in,
                                                             const std::vector<arma::mat> &covariances_in,
                                                             GradientEstimator *gradient_estimator_in)
      : ProposalKernel()
  {
    this->gradient_estimator.reset(gradient_estimator_in);
    this->gradient_estimator->set_proposal(this);
    this->index.reset();

    // make a vector of zero means and a vector of identity matrices

    std::vector<arma::colvec> zero_means;
    std::vector<arma::mat> identity_matrices;
    zero_means.reserve(variable_names_in.size());
    identity_matrices.reserve(variable_names_in.size());

    for (size_t i = 0;
         i < variable_names_in.size();
         ++i)
    {
      zero_means.push_back(arma::zeros<arma::colvec>(covariances_in[i].n_rows));
      identity_matrices.push_back(arma::eye<arma::mat>(covariances_in[i].n_rows, covariances_in[i].n_rows));
      this->proposal_info[variable_names_in[i]] = GaussianProposalInfo(covariances_in[i]);
    }
    this->proposal_simulate.reset(new GaussianIndependentProposalKernel(variable_names_in,
                                                                    zero_means,
                                                                    identity_matrices));
  }

  BarkerDynamicsProposalKernel::BarkerDynamicsProposalKernel(const std::string &variable_name_in,
                                                             const double &sd_in,
                                                             GradientEstimator *gradient_estimator_in)
      : ProposalKernel()
  {
    this->gradient_estimator.reset(gradient_estimator_in);
    this->gradient_estimator->set_proposal(this);
    this->index.reset();

    this->proposal_simulate.reset(new GaussianIndependentProposalKernel(variable_name_in,
                                                                    0.0,
                                                                    1.0));

    this->proposal_info[variable_name_in] = GaussianProposalInfo(sd_in);
  }

  BarkerDynamicsProposalKernel::BarkerDynamicsProposalKernel(const BarkerDynamicsProposalKernel &another)
      : ProposalKernel(another)
  {
    this->make_copy(another);
  }

  void BarkerDynamicsProposalKernel::operator=(const BarkerDynamicsProposalKernel &another)
  {
    if (this == &another)
      return;

    this->gradient_estimator.reset();

    this->proposal_simulate.reset();

    this->index.reset();

    ProposalKernel::operator=(another);
    this->make_copy(another);
  }

  Kernel *BarkerDynamicsProposalKernel::duplicate() const
  {
    return (new BarkerDynamicsProposalKernel(*this));
  }

  ProposalKernel *BarkerDynamicsProposalKernel::proposal_kernel_duplicate() const
  {
    return (new BarkerDynamicsProposalKernel(*this));
  }

  void BarkerDynamicsProposalKernel::make_copy(const BarkerDynamicsProposalKernel &another)
  {
    this->proposal_info = another.proposal_info;
    this->proposal_store = another.proposal_store;
    if (another.gradient_estimator != NULL)
      this->gradient_estimator.reset(another.gradient_estimator->duplicate());
    else
      this->gradient_estimator.reset();
    if (another.proposal_simulate != NULL)
      this->proposal_simulate.reset(another.proposal_simulate->independent_proposal_kernel_duplicate());
    else
      this->proposal_simulate.reset();

    if (another.index != NULL)
      this->index.reset(another.index->duplicate());
    else
      this->index.reset();
  }

  double BarkerDynamicsProposalKernel::specific_evaluate_kernel(const Particle &proposed_particle,
                                                                const Particle &old_particle) const
  {
    GradientEstimatorOutput *gradient_estimator_output = old_particle.get_gradient_estimator_output(this);

    double output = 0.0;
    for (auto i = this->proposal_info.begin();
         i != this->proposal_info.end();
         ++i)
    {
      arma::colvec current_chol_gradient_prod = sqrt(i->second.get_double_scale()) * i->second.get_chol() * arma::vectorise(gradient_estimator_output->get_gradient_of_log(i->first, this->index.get(), old_particle));

      arma::colvec current_proposed = proposed_particle.get_transformed_parameters(this).get_colvec(i->first);
      arma::colvec old_mat = old_particle.get_transformed_parameters(this).get_colvec(i->first);
      // z = (C^T)^{-1}(y - x) with C = chol(s Sigma) = sqrt(s) chol(Sigma)
      // (Livingstone and Zanella, supplement, preconditioned Barker).
      arma::rowvec z = (current_proposed - arma::vectorise(old_mat)).t() * i->second.get_inv_chol() / sqrt(i->second.get_double_scale());

      for (size_t j = 0; j < current_proposed.size(); ++j)
      {
        output = output - log1p(exp(-z[j] * current_chol_gradient_prod[j]));
      }
    }
    return output;
  }

  double BarkerDynamicsProposalKernel::specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                                          const Particle &old_particle) const
  {
    GradientEstimatorOutput *gradient_estimator_output = old_particle.get_gradient_estimator_output(this);

    double output = 0.0;
    for (auto i = this->proposal_info.begin();
         i != this->proposal_info.end();
         ++i)
    {
      arma::colvec current_chol_gradient_prod = sqrt(i->second.get_double_scale()) * i->second.get_chol() * arma::vectorise(gradient_estimator_output->subsample_get_gradient_of_log(i->first, this->index.get(), old_particle));

      arma::colvec current_proposed = proposed_particle.get_transformed_parameters(this).get_colvec(i->first);
      arma::colvec old_mat = old_particle.get_transformed_parameters(this).get_colvec(i->first);
      // z = (C^T)^{-1}(y - x) with C = chol(s Sigma) = sqrt(s) chol(Sigma)
      // (Livingstone and Zanella, supplement, preconditioned Barker).
      arma::rowvec z = (current_proposed - arma::vectorise(old_mat)).t() * i->second.get_inv_chol() / sqrt(i->second.get_double_scale());

      for (size_t j = 0; j < current_proposed.size(); ++j)
      {
        output = output - log1p(exp(-z[j] * current_chol_gradient_prod[j]));
      }
    }
    return output;
  }

  Parameters BarkerDynamicsProposalKernel::simulate(RandomNumberGenerator &rng,
                                                    const Particle &current_particle) const
  {
    GradientEstimatorOutput *gradient_estimator_output = current_particle.get_gradient_estimator_output(this);

    // points simulated are in transformed space (the standard normal steps,
    // under each block's name)
    Parameters proposed = this->proposal_simulate->independent_simulate(rng);
    Parameters output;

    for (auto i = this->proposal_info.begin();
         i != this->proposal_info.end();
         ++i)
    {
      arma::colvec current_chol_gradient_prod = sqrt(i->second.get_double_scale()) * i->second.get_chol() * arma::vectorise(gradient_estimator_output->get_gradient_of_log(i->first, this->index.get(), current_particle));

      arma::rowvec current_proposed = arma::vectorise(proposed[i->first]).t();
      check_proposal_size(i->first, current_particle.get_transformed_parameters(this).get_colvec(i->first).n_elem,
                          i->second.get_chol().n_rows, "covariance");

      arma::rowvec log_unifs = log(multiple_runif(rng, current_proposed.size()));
      for (size_t j = 0; j < current_proposed.size(); ++j)
      {
        double log_prob = -log1p(exp(-current_proposed[j] * current_chol_gradient_prod[j]));
        if (log_unifs[j] > log_prob)
        {
          current_proposed[j] = -current_proposed[j];
        }
      }
      const Parameters &current = current_particle.get_transformed_parameters(this);
      output.set_stacked(i->first, current.get_colvec(i->first) + (current_proposed * sqrt(i->second.get_double_scale()) * i->second.get_chol()).t(), current);
    }

    return output;
  }

  Parameters BarkerDynamicsProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                              const Particle &particle) const
  {
    Rcpp::stop("BarkerDynamicsProposalKernel::subsample_simulate - not written yet.");
  }

  Parameters BarkerDynamicsProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                              const std::string &variable,
                                                              const Particle &proposed_particle) const
  {
    Rcpp::stop("BarkerDynamicsProposalKernel::subsample_simulate - not written yet.");
  }

  arma::mat BarkerDynamicsProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                                   const Particle &proposed_particle,
                                                                   const Particle &old_particle)
  {
    Rcpp::stop("BarkerDynamicsProposalKernel::specific_gradient_of_log - not written yet.");
  }

  arma::mat BarkerDynamicsProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                             const Particle &proposed_particle,
                                                                             const Particle &old_particle)
  {
    Rcpp::stop("BarkerDynamicsProposalKernel::specific_subsample_gradient_of_log - not written yet.");
  }

  void BarkerDynamicsProposalKernel::set_proposal_parameters(Parameters *proposal_parameters_in)
  {
  }

  GradientEstimatorOutput *BarkerDynamicsProposalKernel::simulate_gradient_estimator_output() const
  {
    GradientEstimatorOutput *current_output = gradient_estimator->initialise();
    current_output->simulate_auxiliary_variables();
    return current_output;
  }

  std::vector<const ProposalKernel *> BarkerDynamicsProposalKernel::get_proposals() const
  {
    std::vector<const ProposalKernel *> proposals = this->proposal_simulate->get_proposals();
    proposals.push_back(this);
    return proposals;
  }

  // The kernel keeps its own private copy of the index: the pointer passed in
  // is owned (and deleted) by the calling MCMC object, so storing it directly
  // would lead to it being deleted twice.
  void BarkerDynamicsProposalKernel::set_index(Index *index_in)
  {
    this->index.reset();
    if (index_in != NULL)
      this->index.reset(index_in->duplicate());
    else
      this->index.reset();
  }

  void BarkerDynamicsProposalKernel::set_index_if_null(Index *index_in)
  {
    if ((this->index == NULL) && (index_in != NULL))
      this->index.reset(index_in->duplicate());
  }

  bool BarkerDynamicsProposalKernel::can_be_evaluated() const
  {
    return true;
  }

  void BarkerDynamicsProposalKernel::set_data(Data *data_in)
  {
  }

GaussianProposalInfo* BarkerDynamicsProposalKernel::get_gaussian_proposal_info(const std::string &variable)
{
  auto found = this->proposal_info.find(variable);
  if (found==this->proposal_info.end())
    return NULL;
  return &found->second;
}

}
