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
    this->gradient_estimator = NULL;
    this->proposal_simulate = NULL;
    this->index = NULL;
  }

  BarkerDynamicsProposalKernel::~BarkerDynamicsProposalKernel()
  {
    if (this->gradient_estimator != NULL)
      delete this->gradient_estimator;

    if (this->proposal_simulate != NULL)
      delete this->proposal_simulate;

    if (this->index != NULL)
      delete index;
  }

  /*
    BarkerDynamicsProposalKernel::BarkerDynamicsProposalKernel(const std::vector<std::string> &variable_names_in,
                                                               GradientEstimator *gradient_estimator_in)
        : ProposalKernel()
    {
      this->gradient_estimator = gradient_estimator_in;
      this->gradient_estimator->set_proposal(this);
      this->proposal_simulate = NULL;
      this->index = NULL;

      for (auto i = variable_names_in.begin();
           i != variable_names_in.end();
           ++i)
      {
        this->proposal_info[*i] = GaussianProposalInfo();
      }
    }
  */

  BarkerDynamicsProposalKernel::BarkerDynamicsProposalKernel(const std::string &variable_name_in,
                                                             const arma::mat &covariance_in,
                                                             GradientEstimator *gradient_estimator_in)
      : ProposalKernel()
  {
    this->gradient_estimator = gradient_estimator_in;
    this->gradient_estimator->set_proposal(this);
    this->index = NULL;
    arma::colvec mean = arma::zeros<arma::colvec>(covariance_in.n_rows);
    // Create an identity matrix for the covariance
    arma::mat identity = arma::eye<arma::mat>(covariance_in.n_rows, covariance_in.n_rows);
    this->proposal_simulate = new GaussianIndependentProposalKernel(variable_name_in,
                                                                    mean,
                                                                    identity);

    this->proposal_info[variable_name_in] = GaussianProposalInfo(covariance_in);
  }

  BarkerDynamicsProposalKernel::BarkerDynamicsProposalKernel(const std::string &variable_name_in,
                                                             const arma::mat &covariance_in,
                                                             double scale_in,
                                                             GradientEstimator *gradient_estimator_in)
      : ProposalKernel()
  {
    this->gradient_estimator = gradient_estimator_in;
    this->gradient_estimator->set_proposal(this);
    this->index = NULL;
    arma::colvec mean = arma::zeros<arma::colvec>(covariance_in.n_rows);
    // Create an identity matrix for the covariance
    arma::mat identity = arma::eye<arma::mat>(covariance_in.n_rows, covariance_in.n_rows);
    this->proposal_simulate = new GaussianIndependentProposalKernel(variable_name_in,
                                                                    mean,
                                                                    identity);

    this->proposal_info[variable_name_in] = GaussianProposalInfo(covariance_in,
                                                                 scale_in);
  }

  BarkerDynamicsProposalKernel::BarkerDynamicsProposalKernel(const std::vector<std::string> &variable_names_in,
                                                             const std::vector<arma::mat> &covariances_in,
                                                             GradientEstimator *gradient_estimator_in)
      : ProposalKernel()
  {
    this->gradient_estimator = gradient_estimator_in;
    this->gradient_estimator->set_proposal(this);
    this->index = NULL;

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
    this->proposal_simulate = new GaussianIndependentProposalKernel(variable_names_in,
                                                                    zero_means,
                                                                    identity_matrices);
  }

  BarkerDynamicsProposalKernel::BarkerDynamicsProposalKernel(const std::string &variable_name_in,
                                                             const double &sd_in,
                                                             GradientEstimator *gradient_estimator_in)
      : ProposalKernel()
  {
    this->gradient_estimator = gradient_estimator_in;
    this->gradient_estimator->set_proposal(this);
    this->index = NULL;

    this->proposal_simulate = new GaussianIndependentProposalKernel(variable_name_in,
                                                                    0.0,
                                                                    1.0);

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

    if (this->gradient_estimator != NULL)
      delete this->gradient_estimator;

    if (this->proposal_simulate != NULL)
      delete this->proposal_simulate;

    if (this->index != NULL)
      delete index;

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
      this->gradient_estimator = another.gradient_estimator->duplicate();
    else
      this->gradient_estimator = NULL;
    if (another.proposal_simulate != NULL)
      this->proposal_simulate = another.proposal_simulate->independent_proposal_kernel_duplicate();
    else
      this->proposal_simulate = NULL;

    if (another.index != NULL)
      this->index = another.index->duplicate();
    else
      this->index = NULL;
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
      arma::colvec current_chol_gradient_prod = sqrt(i->second.get_double_scale()) * i->second.get_chol() * arma::vectorise(gradient_estimator_output->get_gradient_of_log(i->first, this->index, old_particle));

      arma::mat proposed_mat = proposed_particle.get_transformed_parameters(this)[i->first];
      arma::mat old_mat = old_particle.get_transformed_parameters(this)[i->first];
      arma::colvec current_proposed = arma::vectorise(proposed_mat);
      arma::rowvec z = (current_proposed - arma::vectorise(old_mat)) * i->second.get_inv_chol();

      for (size_t j = 0; j < current_proposed.size(); ++j)
      {
        output = output - log1p(exp(-z[j] * current_chol_gradient_prod[j]));
      }
    }
    return output;
  }

  /*
   double BarkerDynamicsProposalKernel::specific_evaluate_kernel(Particle &proposed_particle,
   Particle &old_particle,
   const Parameters &conditioned_on_parameters) const
   {
   GradientEstimatorOutput* estimator = old_particle.initialise_gradient_estimator_output(this,
   this->gradient_estimator);

   double output = 0.0;
   for (auto i=this->proposal_info.begin();
   i!=this->proposal_info.end();
   ++i)
   {
   arma::colvec current_chol_gradient_prod = sqrt(i->second.get_double_scale())*i->second.get_chol() * arma::vectorise(estimator->get_gradient_of_log(i->first,this->index,old_particle,conditioned_on_parameters));

   arma::mat proposed_mat = proposed_particle.parameters[i->first];
   arma::mat old_mat = old_particle.parameters[i->first];
   arma::colvec current_proposed = arma::vectorise(proposed_mat);
   arma::rowvec z = (current_proposed - arma::vectorise(old_mat))*i->second.get_inv_chol();

   for (size_t j=0; j<current_proposed.size(); ++j)
   {
   output = output - log1p(exp(-z[j]*current_chol_gradient_prod[j]));
   }
   }
   return output;
   }
   */

  double BarkerDynamicsProposalKernel::specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                                          const Particle &old_particle) const
  {
    GradientEstimatorOutput *gradient_estimator_output = old_particle.get_gradient_estimator_output(this);

    double output = 0.0;
    for (auto i = this->proposal_info.begin();
         i != this->proposal_info.end();
         ++i)
    {
      arma::colvec current_chol_gradient_prod = sqrt(i->second.get_double_scale()) * i->second.get_chol() * arma::vectorise(gradient_estimator_output->subsample_get_gradient_of_log(i->first, this->index, old_particle));

      arma::mat proposed_mat = proposed_particle.get_transformed_parameters(this)[i->first];
      arma::mat old_mat = old_particle.get_transformed_parameters(this)[i->first];
      arma::colvec current_proposed = arma::vectorise(proposed_mat);
      arma::rowvec z = (current_proposed - arma::vectorise(old_mat)) * i->second.get_inv_chol();

      for (size_t j = 0; j < current_proposed.size(); ++j)
      {
        output = output - log1p(exp(-z[j] * current_chol_gradient_prod[j]));
      }
    }
    return output;
  }

  /*
   double BarkerDynamicsProposalKernel::specific_subsample_evaluate_kernel(Particle &proposed_particle,
   Particle &old_particle,
   const Parameters &conditioned_on_parameters) const
   {
   GradientEstimatorOutput* estimator = old_particle.initialise_gradient_estimator_output(this,
   this->gradient_estimator);

   double output = 0.0;
   for (auto i=this->proposal_info.begin();
   i!=this->proposal_info.end();
   ++i)
   {
   arma::colvec current_chol_gradient_prod = sqrt(i->second.get_double_scale())*i->second.get_chol() * arma::vectorise(estimator->subsample_get_gradient_of_log(i->first,this->index,old_particle));

   arma::mat proposed_mat = proposed_particle.parameters[i->first];
   arma::mat old_mat = old_particle.parameters[i->first];
   arma::colvec current_proposed = arma::vectorise(proposed_mat);
   arma::rowvec z = (current_proposed - arma::vectorise(old_mat))*i->second.get_inv_chol();

   for (size_t j=0; j<current_proposed.size(); ++j)
   {
   output = output - log1p(exp(-z[j]*current_chol_gradient_prod[j]));
   }
   }
   return output;
   }
   */

  Parameters BarkerDynamicsProposalKernel::simulate(RandomNumberGenerator &rng,
                                                    const Particle &current_particle) const
  {
    GradientEstimatorOutput *gradient_estimator_output = current_particle.get_gradient_estimator_output(this);

    // Parameters previous_parameters = proposed_particle.previous_proposal_store.find(this)->get_transformed_parameters();

    // points simulated are in transformed space
    Parameters proposed = this->proposal_simulate->simulate(rng, current_particle);

    for (auto i = this->proposal_info.begin();
         i != this->proposal_info.end();
         ++i)
    {
      arma::colvec current_chol_gradient_prod = sqrt(i->second.get_double_scale()) * i->second.get_chol() * arma::vectorise(gradient_estimator_output->get_gradient_of_log(i->first, this->index, current_particle));

      arma::mat initial_proposed = proposed[i->first];
      arma::rowvec current_proposed = arma::vectorise(initial_proposed);

      arma::rowvec log_unifs = log(multiple_runif(rng, current_proposed.size()));
      for (size_t j = 0; j < current_proposed.size(); ++j)
      {
        double log_prob = -log1p(exp(-current_proposed[j] * current_chol_gradient_prod[j]));
        if (log_unifs[j] > log_prob)
        {
          current_proposed[j] = -current_proposed[j];
        }
      }
      proposed[i->first] = current_particle.get_transformed_parameters(this)[i->first] + arma::reshape(current_proposed * sqrt(i->second.get_double_scale()) * i->second.get_chol(), initial_proposed.n_rows, initial_proposed.n_cols);
    }

    return proposed;
  }

  /*
   Parameters BarkerDynamicsProposalKernel::simulate(RandomNumberGenerator &rng,
   Particle &particle,
   const Parameters &conditioned_on_parameters) const
   {
   GradientEstimatorOutput* estimator = particle.initialise_gradient_estimator_output(this,
   this->gradient_estimator);

   Parameters proposed = this->proposal_simulate->simulate(rng,
   particle,
   conditioned_on_parameters);

   for (auto i=this->proposal_info.begin();
   i!=this->proposal_info.end();
   ++i)
   {
   arma::colvec current_chol_gradient_prod = sqrt(i->second.get_double_scale())*i->second.get_chol() * arma::vectorise(estimator->get_gradient_of_log(i->first,this->index,particle,conditioned_on_parameters));

   arma::mat initial_proposed = proposed[i->first];
   arma::rowvec current_proposed = arma::vectorise(initial_proposed);

   arma::rowvec log_unifs = log(multiple_runif(rng,current_proposed.size()));
   for (size_t j=0; j<current_proposed.size(); ++j)
   {
   double log_prob = -log1p(exp(-current_proposed[j]*current_chol_gradient_prod[j]));
   if (log_unifs[j]>log_prob)
   {
   current_proposed[j] = -current_proposed[j];
   }
   }
   proposed[i->first] = particle.parameters[i->first] + arma::reshape(current_proposed * sqrt(i->second.get_double_scale())*i->second.get_chol(),initial_proposed.n_rows,initial_proposed.n_cols);
   }
   return proposed;
   }
   */

  Parameters BarkerDynamicsProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                              const Particle &particle) const
  {
    Rcpp::stop("BarkerDynamicsProposalKernel::subsample_simulate - not written yet.");
    /*
     Parameters proposed = this->proposal_simulate->subsample_simulate(rng,
     particle);

     for (auto i=this->proposal_info.begin();
     i!=this->proposal_info.end();
     ++i)
     {
     arma::colvec current_chol_gradient_prod = i->second.chol * arma::vectorise(this->subsample_gradient_estimator->get_gradient_of_log(i->first,particle));

     arma::mat initial_proposed = proposed[i->first];
     arma::rowvec current_proposed = arma::vectorise(initial_proposed);

     arma::rowvec log_unifs = log(runif(rng,current_proposed.size()));
     for (size_t j=0; j<current_proposed.size(); ++j)
     {
     double log_prob = -log1p(exp(-current_proposed[j]*current_chol_gradient_prod[j]));
     if (log_unifs[j]>log_prob)
     {
     current_proposed[j] = -current_proposed[j];
     }
     }
     proposed[i->first] = particle.parameters[i->first] + arma::reshape(current_proposed * i->second.chol,initial_proposed.n_rows,initial_proposed.n_cols);
     }
     return proposed;
     */
  }

  /*
   Parameters BarkerDynamicsProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
   Particle &particle,
   const Parameters &conditioned_on_parameters) const
   {
   Rcpp::stop("BarkerDynamicsProposalKernel::subsample_simulate - not written yet.");
   }
   */

  Parameters BarkerDynamicsProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                              const std::string &variable,
                                                              const Particle &proposed_particle) const
  {
    Rcpp::stop("BarkerDynamicsProposalKernel::subsample_simulate - not written yet.");
    /*
     Parameters proposed = this->proposal_simulate->subsample_simulate(rng,
     variable,
     particle);

     auto found = this->proposal_info.find(variable);

     arma::colvec current_chol_gradient_prod = found->second.chol * arma::vectorise(this->subsample_gradient_estimator->get_gradient_of_log(found->first,particle));

     arma::mat initial_proposed = proposed[found->first];
     arma::rowvec current_proposed = arma::vectorise(initial_proposed);

     arma::rowvec log_unifs = log(runif(rng,current_proposed.size()));
     for (size_t j=0; j<current_proposed.size(); ++j)
     {
     double log_prob = -log1p(exp(-current_proposed[j]*current_chol_gradient_prod[j]));
     if (log_unifs[j]>log_prob)
     {
     current_proposed[j] = -current_proposed[j];
     }
     }
     proposed[found->first] = particle.parameters[found->first] + arma::reshape(current_proposed * found->second.chol,initial_proposed.n_rows,initial_proposed.n_cols);
     return proposed;
     */
  }

  /*
   Parameters BarkerDynamicsProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
   const std::string &variable,
   Particle &particle,
   const Parameters &conditioned_on_parameters) const
   {
   Rcpp::stop("BarkerDynamicsProposalKernel::subsample_simulate - not written yet.");
   }
   */

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

  /*
   arma::mat BarkerDynamicsProposalKernel::specific_gradient_of_log(const std::string &variable,
   Particle &proposed_particle,
   Particle &old_particle,
   const Parameters &conditioned_on_parameters)
   {
   Rcpp::stop("BarkerDynamicsProposalKernel::specific_gradient_of_log - not written yet.");
   }
   */

  // virtual arma::mat specific_subsample_gradient_of_log(const std::string &variable,
  //                                                      Particle &proposed_particle,
  //                                                      Particle &old_particle)=0;

  /*
   arma::mat BarkerDynamicsProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
   Particle &proposed_particle,
   Particle &old_particle,
   const Parameters &conditioned_on_parameters)
   {
   Rcpp::stop("BarkerDynamicsProposalKernel::specific_gradient_of_log - not written yet.");
   }
   */

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

  void BarkerDynamicsProposalKernel::set_index(Index *index_in)
  {
    this->index = index_in;
  }

  void BarkerDynamicsProposalKernel::set_index_if_null(Index *index_in)
  {
    if (this->index == NULL)
      this->index = index_in;
  }

  bool BarkerDynamicsProposalKernel::can_be_evaluated() const
  {
    return true;
  }

  void BarkerDynamicsProposalKernel::set_data(Data *data_in)
  {
  }
}
