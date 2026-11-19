#include "element_adaptors.h"
#include "element_random_walk_proposal_kernel.h"
#include "covariance_tools.h"
#include "smc_output.h"
#include "particles.h"
#include "move_output.h"
#include "ensemble_kalman_output.h"
#include "ensemble.h"

namespace ilike
{

// --- ElementScaleAdaptor ---------------------------------------------------

ElementScaleAdaptor::ElementScaleAdaptor()
:MCMCAdaptor()
{
  this->target = 0.44;
  this->gain_exponent = 0.6;
  this->gain_constant = 1.0;
  this->averaging = false;
  this->number_of_updates = 0;
  this->log_scale = 0.0;
  this->average_log_scale = 0.0;
  this->started = false;
}

ElementScaleAdaptor::ElementScaleAdaptor(double target_in,
                                         double gain_exponent_in,
                                         double gain_constant_in,
                                         bool averaging_in)
:MCMCAdaptor()
{
  this->target = target_in;
  this->gain_exponent = gain_exponent_in;
  this->gain_constant = gain_constant_in;
  this->averaging = averaging_in;
  this->number_of_updates = 0;
  this->log_scale = 0.0;
  this->average_log_scale = 0.0;
  this->started = false;
}

ElementScaleAdaptor::~ElementScaleAdaptor()
{
}

ElementScaleAdaptor::ElementScaleAdaptor(const ElementScaleAdaptor &another)
:MCMCAdaptor(another)
{
  this->make_copy(another);
}

void ElementScaleAdaptor::operator=(const ElementScaleAdaptor &another)
{
  if (this == &another)
    return;
  MCMCAdaptor::operator=(another);
  this->make_copy(another);
}

MCMCAdaptor* ElementScaleAdaptor::duplicate() const
{
  return new ElementScaleAdaptor(*this);
}

void ElementScaleAdaptor::make_copy(const ElementScaleAdaptor &another)
{
  this->target = another.target;
  this->gain_exponent = another.gain_exponent;
  this->gain_constant = another.gain_constant;
  this->averaging = another.averaging;
  this->number_of_updates = another.number_of_updates;
  this->log_scale = another.log_scale;
  this->average_log_scale = another.average_log_scale;
  this->started = another.started;
}

ElementRandomWalkProposalKernel* ElementScaleAdaptor::kernel() const
{
  ElementRandomWalkProposalKernel* element = dynamic_cast<ElementRandomWalkProposalKernel*>(this->proposal);
  if (element == NULL)
    Rcpp::stop("ElementScaleAdaptor: not attached to a move on an element.");
  return element;
}

void ElementScaleAdaptor::adapt(const AdaptationStep &step)
{
  ElementRandomWalkProposalKernel* element = this->kernel();
  if (!this->started)
  {
    this->log_scale = std::log(element->get_scale());
    this->average_log_scale = this->log_scale;
    this->started = true;
  }
  this->number_of_updates = this->number_of_updates + 1;
  double gain = this->gain_constant*std::pow(double(this->number_of_updates + 1), -this->gain_exponent);
  this->log_scale = this->log_scale + gain*(step.acceptance_probability - this->target);
  if (this->averaging)
    this->average_log_scale = this->average_log_scale +
                              (this->log_scale - this->average_log_scale)/double(this->number_of_updates + 1);
  element->set_scale(std::exp(this->averaging ? this->average_log_scale : this->log_scale));
}

// for a stopping criterion (C207): the element's step, and its target
bool ElementScaleAdaptor::watched(const std::string &name, std::vector<double> &values) const
{
  if (name != "scale")
    return false;
  values.assign(1, this->kernel()->get_scale());
  return true;
}

std::vector<std::string> ElementScaleAdaptor::watchable() const
{
  return std::vector<std::string>(1, "scale");
}

bool ElementScaleAdaptor::acceptance_target(double &target) const
{
  target = this->target;
  return true;
}

std::string ElementScaleAdaptor::state_label() const
{
  ElementRandomWalkProposalKernel* element = this->kernel();
  std::string index = std::to_string(element->get_row() + 1);
  if (element->get_column() > 0)
    index = index + "," + std::to_string(element->get_column() + 1);
  return element->get_variable() + "[" + index + "]";
}

void ElementScaleAdaptor::append_state(std::vector<std::string> *names,
                                       std::vector<double> *values,
                                       bool final,
                                       const std::string &prefix) const
{
  ElementRandomWalkProposalKernel* element = this->kernel();
  if (names != NULL)
  {
    std::string index = std::to_string(element->get_row() + 1);
    if (element->get_column() > 0)
      index = index + "," + std::to_string(element->get_column() + 1);
    names->push_back(prefix + element->get_variable() + "[" + index + "].scale");
  }
  if (values != NULL)
    values->push_back(element->get_scale());
}

// --- ElementSMCAdaptor -----------------------------------------------------

ElementSMCAdaptor::ElementSMCAdaptor()
:SMCAdaptor()
{
  this->constant = 2.38;
}

ElementSMCAdaptor::ElementSMCAdaptor(double constant_in)
:SMCAdaptor()
{
  this->constant = constant_in;
}

ElementSMCAdaptor::~ElementSMCAdaptor()
{
}

ElementSMCAdaptor::ElementSMCAdaptor(const ElementSMCAdaptor &another)
:SMCAdaptor(another)
{
  this->make_copy(another);
}

void ElementSMCAdaptor::operator=(const ElementSMCAdaptor &another)
{
  if (this == &another)
    return;
  SMCAdaptor::operator=(another);
  this->make_copy(another);
}

SMCAdaptor* ElementSMCAdaptor::duplicate() const
{
  return new ElementSMCAdaptor(*this);
}

void ElementSMCAdaptor::make_copy(const ElementSMCAdaptor &another)
{
  this->constant = another.constant;
}

void ElementSMCAdaptor::adapt_from(const arma::mat &points,
                                   const arma::colvec &normalised_weights,
                                   const std::vector<std::string> &names,
                                   const std::vector<size_t> &offsets,
                                   size_t number_of_rows)
{
  ElementRandomWalkProposalKernel* element = dynamic_cast<ElementRandomWalkProposalKernel*>(this->proposal);
  if (element == NULL)
    Rcpp::stop("ElementSMCAdaptor: not attached to a move on an element.");
  auto found = std::find(names.begin(), names.end(), element->get_variable());
  if ((found == names.end()) || (points.n_rows < 2))
    return;
  size_t k = found - names.begin();
  arma::uword column = offsets[k] + element->get_column()*number_of_rows + element->get_row();
  if (column >= points.n_cols)
    return;
  // the element first, then every other column
  arma::uvec order(points.n_cols);
  order[0] = column;
  arma::uword next = 1;
  for (arma::uword c = 0; c < points.n_cols; ++c)
    if (c != column)
      order[next++] = c;
  arma::rowvec mean;
  arma::mat covariance;
  if (!conditional_moments(points.cols(order), normalised_weights, 1, mean, covariance))
    return;
  double sd = std::sqrt(covariance(0, 0));
  if (std::isfinite(sd) && (sd > 0.0))
    element->set_scale(this->constant*sd);
}

// the number of rows of the element's variable, at a point
static size_t rows_of(const Parameters &point, const std::string &variable)
{
  return point[variable].n_rows;
}

void ElementSMCAdaptor::smc_adapt(SMCOutput* current_state)
{
  ElementRandomWalkProposalKernel* element = dynamic_cast<ElementRandomWalkProposalKernel*>(this->proposal);
  if (element == NULL)
    Rcpp::stop("ElementSMCAdaptor: not attached to a move on an element.");
  Particles &generation = current_state->back();
  std::vector<std::string> names;
  std::vector<size_t> offsets;
  arma::mat points = population_points(generation, names, offsets);
  if (points.n_rows == 0)
    return;
  size_t number_of_rows = generation.is_flat() ? rows_of(generation.flat_points.parameters_at(0), element->get_variable())
                                               : rows_of(generation.particles[0]->back().parameters, element->get_variable());
  arma::colvec weights;
  if (!population_weights(generation, points.n_rows, weights))
    return;
  this->adapt_from(points, weights, names, offsets, number_of_rows);
}

void ElementSMCAdaptor::ensemble_adapt(EnsembleKalmanOutput* current_state)
{
  ElementRandomWalkProposalKernel* element = dynamic_cast<ElementRandomWalkProposalKernel*>(this->proposal);
  if (element == NULL)
    Rcpp::stop("ElementSMCAdaptor: not attached to a move on an element.");
  Ensemble &ensemble = current_state->back();
  std::vector<std::string> names;
  std::vector<size_t> offsets;
  arma::mat points = ensemble_points(ensemble, names, offsets);
  if (points.n_rows == 0)
    return;
  size_t number_of_rows = rows_of(ensemble.members[0]->back().parameters, element->get_variable());
  this->adapt_from(points, arma::colvec(points.n_rows, arma::fill::value(1.0/double(points.n_rows))),
                   names, offsets, number_of_rows);
}

}
