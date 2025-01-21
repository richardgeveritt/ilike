#include "truncated_gaussian_distribution_factor.h"
#include <iterator>

namespace ilike {
TruncatedGaussianDistributionFactor::TruncatedGaussianDistributionFactor()
    : DistributionFactor() {}

TruncatedGaussianDistributionFactor::TruncatedGaussianDistributionFactor(
    const std::string &variable_in, double mean_in, double sd_in,
    double lower_in, double upper_in)
    : DistributionFactor() {
  this->proposal_info[variable_in] = GaussianProposalInfo(mean_in, sd_in);
  this->lower_info[variable_in] = lower_in;
  this->upper_info[variable_in] = upper_in;
}

TruncatedGaussianDistributionFactor::TruncatedGaussianDistributionFactor(
    const std::string &variable_in, const arma::colvec &mean_in,
    const arma::mat &covariance_in, const arma::colvec &lower_in,
    const arma::colvec &upper_in)
    : DistributionFactor() {
  this->proposal_info[variable_in] =
      GaussianProposalInfo(mean_in, covariance_in);
  this->lower_info[variable_in] = lower_in;
  this->upper_info[variable_in] = upper_in;
}

TruncatedGaussianDistributionFactor::TruncatedGaussianDistributionFactor(
    const std::vector<std::string> &variable_names_in,
    const std::vector<arma::colvec> &means_in,
    const std::vector<arma::mat> &covariances_in,
    const std::vector<arma::colvec> &lower_in,
    const std::vector<arma::colvec> &upper_in)
    : DistributionFactor() {
  for (size_t i = 0; i < variable_names_in.size(); ++i) {
    this->proposal_info[variable_names_in[i]] =
        GaussianProposalInfo(means_in[i], covariances_in[i]);
    this->lower_info[variable_names_in[i]] = lower_in[i];
    this->upper_info[variable_names_in[i]] = upper_in[i];
  }
}

TruncatedGaussianDistributionFactor::~TruncatedGaussianDistributionFactor() {}

TruncatedGaussianDistributionFactor::TruncatedGaussianDistributionFactor(
    const TruncatedGaussianDistributionFactor &another)
    : DistributionFactor(another) {
  this->make_copy(another);
}

void TruncatedGaussianDistributionFactor::operator=(
    const TruncatedGaussianDistributionFactor &another) {
  if (this == &another)
    return;

  DistributionFactor::operator=(another);
  this->make_copy(another);
}

Factor *TruncatedGaussianDistributionFactor::duplicate() const {
  return (new TruncatedGaussianDistributionFactor(*this));
}

DistributionFactor *
TruncatedGaussianDistributionFactor::distribution_factor_duplicate() const {
  return (new TruncatedGaussianDistributionFactor(*this));
}

void TruncatedGaussianDistributionFactor::make_copy(
    const TruncatedGaussianDistributionFactor &another) {
  this->proposal_info = another.proposal_info;
  this->lower_info = another.lower_info;
  this->upper_info = another.upper_info;
}

double TruncatedGaussianDistributionFactor::distribution_evaluate(
    const Parameters &input) const {
      
      Rcpp::stop("TruncatedGaussianDistributionFactor::distribution_evaluate - not written yet.");
      /*
  double output = 0.0;
  for (auto i = this->proposal_info.begin(); i != this->proposal_info.end();
       ++i) {
    arma::colvec mean = i->second.get_mean();
    double scale = i->second.get_double_scale();
    double dim = double(mean.n_rows);
    output = output +
             dtmvnorm_using_precomp(input.get_colvec(i->first), mean,
                                    (1.0 / scale) * i->second.get_inv(),
                                    dim * log(scale) + i->second.get_logdet(),
                                    this->lower_info.at(i->first),
                                    this->upper_info.at(i->first));
  }
  return output;
       */
}

arma::mat TruncatedGaussianDistributionFactor::distribution_evaluate_gradient(
    const std::string &variable, const Parameters &input) const {
  Rcpp::stop("TruncatedGaussianDistributionFactor::distribution_evaluate_"
             "gradient - not written yet.");
  // auto found = this->proposal_info.find(variable);
  // return -(1.0 / found->second.get_double_scale()) * found->second.get_inv()
  // * (input[variable] - found->second.get_mean());
}
} // namespace ilike
