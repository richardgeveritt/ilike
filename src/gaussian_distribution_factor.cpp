#include <iterator>
#include "gaussian_distribution_factor.h"

GaussianDistributionFactor::GaussianDistributionFactor()
  :DistributionFactor()
{
}

GaussianDistributionFactor::GaussianDistributionFactor(const std::string &variable_in,
                                                       double mean_in,
                                                       double sd_in)
:DistributionFactor()
{
    this->proposal_info[variable_in] = GaussianProposalInfo(mean_in,sd_in);
}

GaussianDistributionFactor::GaussianDistributionFactor(const std::string &variable_in,
                                                       const arma::colvec &mean_in,
                                                       const arma::mat &covariance_in)
:DistributionFactor()
{
  this->proposal_info[variable_in] = GaussianProposalInfo(mean_in,covariance_in);
}

GaussianDistributionFactor::GaussianDistributionFactor(const std::vector<std::string> &variable_names_in,
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

GaussianDistributionFactor::~GaussianDistributionFactor()
{
}

GaussianDistributionFactor::GaussianDistributionFactor(const GaussianDistributionFactor &another)
  :DistributionFactor(another)
{
  this->make_copy(another);
}

void GaussianDistributionFactor::operator=(const GaussianDistributionFactor &another)
{
  if(this == &another)
    return;

  DistributionFactor::operator=(another);
  this->make_copy(another);
}

Factor* GaussianDistributionFactor::duplicate() const
{
  return( new GaussianDistributionFactor(*this));
}

DistributionFactor* GaussianDistributionFactor::distribution_factor_duplicate() const
{
  return( new GaussianDistributionFactor(*this));
}

void GaussianDistributionFactor::make_copy(const GaussianDistributionFactor &another)
{
  this->proposal_info = another.proposal_info;
}

double GaussianDistributionFactor::distribution_evaluate(const Parameters &input) const
{
  double output = 0.0;
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::colvec mean = i->second.get_mean();
    double scale = i->second.get_double_scale();
    double dim = double(mean.n_rows);
    output = output + dmvnorm_using_precomp(input.get_colvec(i->first),
                                            mean,
                                            (1.0/scale)*i->second.get_inv(),
                                            dim*log(scale)+i->second.get_logdet());
  }
  return output;
}

arma::mat GaussianDistributionFactor::distribution_evaluate_gradient(const std::string &variable,
                                                                     const Parameters &input) const
{
  auto found = this->proposal_info.find(variable);
  return -(1.0/found->second.get_double_scale())*found->second.get_inv()*(input[variable]-found->second.get_mean());
}
