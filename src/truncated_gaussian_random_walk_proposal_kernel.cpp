#include "truncated_gaussian_random_walk_proposal_kernel.h"
#include "distributions.h"
#include "likelihood_estimator.h"
#include "likelihood_estimator_output.h"
#include <iterator>

namespace ilike {
TruncatedGaussianRandomWalkProposalKernel::
    TruncatedGaussianRandomWalkProposalKernel()
    : SymmetricProposalKernel() {
}

TruncatedGaussianRandomWalkProposalKernel::
    ~TruncatedGaussianRandomWalkProposalKernel() {}

TruncatedGaussianRandomWalkProposalKernel::
    TruncatedGaussianRandomWalkProposalKernel(
        const std::vector<std::string> &variable_names_in,
        const std::vector<arma::colvec> &lower_in,
        const std::vector<arma::colvec> &upper_in)
    : SymmetricProposalKernel() {
  size_t counter = 0;
  for (auto i = variable_names_in.begin(); i != variable_names_in.end(); ++i) {
    this->proposal_info[*i] = GaussianProposalInfo();
    this->lower_info[*i] = lower_in[counter];
    this->upper_info[*i] = upper_in[counter];
    counter = counter + 1;
  }
}

TruncatedGaussianRandomWalkProposalKernel::
    TruncatedGaussianRandomWalkProposalKernel(
        const std::vector<std::string> &variable_names_in,
        const std::vector<arma::mat> &covariances_in,
        const std::vector<arma::colvec> &lower_in,
        const std::vector<arma::colvec> &upper_in)
    : SymmetricProposalKernel() {
  for (size_t i = 0; i < variable_names_in.size(); ++i) {
    this->proposal_info[variable_names_in[i]] =
        GaussianProposalInfo(covariances_in[i]);
    this->lower_info[variable_names_in[i]] = lower_in[i];
    this->upper_info[variable_names_in[i]] = upper_in[i];
  }
}

TruncatedGaussianRandomWalkProposalKernel::
    TruncatedGaussianRandomWalkProposalKernel(
        const std::string &variable_name_in, const double &sd_in,
        const double &lower_in, const double &upper_in)
    : SymmetricProposalKernel() {
  this->proposal_info[variable_name_in] = GaussianProposalInfo(sd_in);
  this->lower_info[variable_name_in] = lower_in;
  this->upper_info[variable_name_in] = upper_in;
}

TruncatedGaussianRandomWalkProposalKernel::
    TruncatedGaussianRandomWalkProposalKernel(
        const std::string &variable_name_in, const arma::mat &covariance_in,
        const arma::colvec &lower_in, const arma::colvec &upper_in)
    : SymmetricProposalKernel() {
  this->proposal_info[variable_name_in] = GaussianProposalInfo(covariance_in);
  this->lower_info[variable_name_in] = lower_in;
  this->upper_info[variable_name_in] = upper_in;
}

TruncatedGaussianRandomWalkProposalKernel::
    TruncatedGaussianRandomWalkProposalKernel(
        const std::string &variable_name_in, const arma::mat &covariance_in,
        double scale_in, const arma::colvec &lower_in,
        const arma::colvec &upper_in)
    : SymmetricProposalKernel() {
  this->proposal_info[variable_name_in] =
      GaussianProposalInfo(covariance_in, scale_in);
  this->lower_info[variable_name_in] = lower_in;
  this->upper_info[variable_name_in] = upper_in;
}

TruncatedGaussianRandomWalkProposalKernel::
    TruncatedGaussianRandomWalkProposalKernel(
        const TruncatedGaussianRandomWalkProposalKernel &another)
    : SymmetricProposalKernel(another) {
  this->make_copy(another);
}

void TruncatedGaussianRandomWalkProposalKernel::operator=(
    const TruncatedGaussianRandomWalkProposalKernel &another) {
  if (this == &another)
    return;

  SymmetricProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel *TruncatedGaussianRandomWalkProposalKernel::duplicate() const {
  return (new TruncatedGaussianRandomWalkProposalKernel(*this));
}

ProposalKernel *
TruncatedGaussianRandomWalkProposalKernel::proposal_kernel_duplicate() const {
  return (new TruncatedGaussianRandomWalkProposalKernel(*this));
}

SymmetricProposalKernel *
TruncatedGaussianRandomWalkProposalKernel::symmetric_proposal_kernel_duplicate()
    const {
  return (new TruncatedGaussianRandomWalkProposalKernel(*this));
}

void TruncatedGaussianRandomWalkProposalKernel::make_copy(
    const TruncatedGaussianRandomWalkProposalKernel &another) {
  this->proposal_info = another.proposal_info;
  this->lower_info = another.lower_info;
  this->upper_info = another.upper_info;
}

double TruncatedGaussianRandomWalkProposalKernel::specific_evaluate_kernel(
    const Particle &proposed_particle, const Particle &old_particle) const {
      Rcpp::stop("TruncatedGaussianDistributionFactor::distribution_evaluate - not written yet.");
      
}

double
TruncatedGaussianRandomWalkProposalKernel::specific_subsample_evaluate_kernel(
    const Particle &proposed_particle, const Particle &old_particle) const {
  // no difference since size of data set does not impact on proposal
  return this->specific_evaluate_kernel(proposed_particle, old_particle);
}

void TruncatedGaussianRandomWalkProposalKernel::set_covariance(
    const std::string &variable, const arma::mat &covariance_in) {
  this->proposal_info[variable].set_covariance(covariance_in);
}

arma::mat TruncatedGaussianRandomWalkProposalKernel::get_inverse_covariance(
    const std::string &variable) {
  return this->proposal_info[variable].get_inv();
}

arma::mat TruncatedGaussianRandomWalkProposalKernel::get_covariance(
    const std::string &variable) {
  return this->proposal_info[variable].get_covariance();
}

Parameters TruncatedGaussianRandomWalkProposalKernel::simulate(
    RandomNumberGenerator &rng, const Particle &particle) const {
  Parameters output;
  for (auto i = this->proposal_info.begin(); i != this->proposal_info.end();
       ++i) {
    double scale = i->second.get_double_scale();
    output[i->first] = rtmvnorm_using_chol(
        rng, particle.get_transformed_parameters(this).get_colvec(i->first),
        sqrt(scale) * i->second.get_chol(), this->lower_info.at(i->first),
        this->upper_info.at(i->first));
  }
  return output;
}

Parameters TruncatedGaussianRandomWalkProposalKernel::subsample_simulate(
    RandomNumberGenerator &rng, const Particle &particle) const {
  // no difference since size of data set does not impact on proposal
  return this->simulate(rng, particle);
}

Parameters TruncatedGaussianRandomWalkProposalKernel::subsample_simulate(
    RandomNumberGenerator &rng, const std::string &variable,
    const Particle &particle) const {
  // no difference since size of data set does not impact on proposal
  auto found = this->proposal_info.find(variable);

  Parameters output;
  output[variable] = rtmvnorm_using_chol(
      rng, particle.get_transformed_parameters(this).get_colvec(variable),
      sqrt(found->second.get_double_scale()) * found->second.get_chol(),
      this->lower_info.at(variable), this->upper_info.at(variable));
  return output;
}

arma::mat TruncatedGaussianRandomWalkProposalKernel::specific_gradient_of_log(
    const std::string &variable, const Particle &proposed_particle,
    const Particle &old_particle) {
  Rcpp::stop("TruncatedGaussianRandomWalkProposalKernel::specific_gradient_of_"
             "log - not written yet.");
}

arma::mat
TruncatedGaussianRandomWalkProposalKernel::specific_subsample_gradient_of_log(
    const std::string &variable, const Particle &proposed_particle,
    const Particle &old_particle) {
  Rcpp::stop("TruncatedGaussianRandomWalkProposalKernel::specific_gradient_of_"
             "log - not written yet.");
}

void TruncatedGaussianRandomWalkProposalKernel::set_proposal_parameters(
    Parameters *proposal_parameters_in) {}

std::vector<std::string>
TruncatedGaussianRandomWalkProposalKernel::get_variables() const {
  std::vector<std::string> variables;
  variables.reserve(this->proposal_info.size());
  for (auto i = this->proposal_info.begin(); i != this->proposal_info.end();
       ++i) {
    variables.push_back(i->first);
  }
  return variables;
}

GradientEstimatorOutput *
TruncatedGaussianRandomWalkProposalKernel::simulate_gradient_estimator_output()
    const {
  return NULL;
}

std::vector<const ProposalKernel *>
TruncatedGaussianRandomWalkProposalKernel::get_proposals() const {
  std::vector<const ProposalKernel *> output;
  output.push_back(this);
  return output;
}

void TruncatedGaussianRandomWalkProposalKernel::set_index(Index *index_in) {}

void TruncatedGaussianRandomWalkProposalKernel::set_index_if_null(
    Index *index_in) {}

bool TruncatedGaussianRandomWalkProposalKernel::can_be_evaluated() const {
  return true;
}

void TruncatedGaussianRandomWalkProposalKernel::set_data(Data *data_in) {}
} // namespace ilike
