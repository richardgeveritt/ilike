#include <iterator>
#include "loggaussian_distribution_factor.h"

LogGaussianDistributionFactor::LogGaussianDistributionFactor()
:DistributionFactor()
{
}

LogGaussianDistributionFactor::LogGaussianDistributionFactor(const std::string &variable_name_in,
                                                             const double &mean_in,
                                                             const double &covariance_in)
:DistributionFactor()
{
  this->proposal_info[variable_name_in] = GaussianProposalInfo(mean_in,covariance_in);
}

LogGaussianDistributionFactor::LogGaussianDistributionFactor(const std::string &variable_name_in,
                                                             const arma::colvec &mean_in,
                                                             const arma::mat &covariance_in)
:DistributionFactor()
{
  this->proposal_info[variable_name_in] = GaussianProposalInfo(mean_in,covariance_in);
}

LogGaussianDistributionFactor::LogGaussianDistributionFactor(const std::vector<std::string> &variable_names_in,
                                                             const std::vector<arma::colvec> &means_in,
                                                             const std::vector<arma::mat> &covariances_in)
:DistributionFactor()
{
  for (size_t i=0;
       i<variable_names_in.size();
       ++i)
  {
    this->proposal_info[variable_names_in[i]] = GaussianProposalInfo(means_in[i],covariances_in[i]);
  }
}

LogGaussianDistributionFactor::~LogGaussianDistributionFactor()
{
}

LogGaussianDistributionFactor::LogGaussianDistributionFactor(const LogGaussianDistributionFactor &another)
:DistributionFactor(another)
{
  this->make_copy(another);
}

void LogGaussianDistributionFactor::operator=(const LogGaussianDistributionFactor &another)
{
  if(this == &another)
    return;
  
  DistributionFactor::operator=(another);
  this->make_copy(another);
}

Factor* LogGaussianDistributionFactor::duplicate() const
{
  return( new LogGaussianDistributionFactor(*this));
}

DistributionFactor* LogGaussianDistributionFactor::distribution_factor_duplicate() const
{
  return( new LogGaussianDistributionFactor(*this));
}

void LogGaussianDistributionFactor::make_copy(const LogGaussianDistributionFactor &another)
{
  this->proposal_info = another.proposal_info;
}

double LogGaussianDistributionFactor::distribution_evaluate(const Parameters &input) const
{
  double output = 0.0;
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::colvec mean = i->second.get_mean();
    double scale = i->second.get_double_scale();
    double dim = double(mean.n_rows);
    output = output + dmvlnorm_using_precomp(input.get_colvec(i->first),
                                             mean,
                                             (1.0/scale)*i->second.get_inv(),
                                             dim*log(scale)+i->second.get_logdet());
  }
  return output;
}

arma::mat LogGaussianDistributionFactor::distribution_evaluate_gradient(const std::string &variable,
                                                                        const Parameters &input) const
{
  //auto found = this->proposal_info.find(variable);
  Rcpp::stop("LogGaussianDistributionFactor::distribution_evaluate_gradient - not written yet.");
}
