#include "plate_likelihood_estimator.h"
#include "plate_likelihood_estimator_output.h"

namespace ilike
{

PlateLikelihoodEstimator::PlateLikelihoodEstimator()
  :ExactLikelihoodEstimator()
{
  this->structure = NULL;
  this->term_function = NULL;
  this->check = false;
  this->current_plate_data = NULL;
  this->structure_built = false;
  this->structure_version = 0;
}

PlateLikelihoodEstimator::PlateLikelihoodEstimator(RandomNumberGenerator* rng_in,
                                                   size_t* seed_in,
                                                   Data* data_in,
                                                   PlateStructurePtr structure_in,
                                                   PlateTermPtr term_in,
                                                   const std::vector<std::string> &variables_in,
                                                   const std::string &name_in,
                                                   bool check_in,
                                                   bool smcfixed_flag_in)
  :ExactLikelihoodEstimator(rng_in, seed_in, data_in,
                            std::vector<DistributionFactor*>(),
                            std::vector<LikelihoodFactor*>(),
                            smcfixed_flag_in)
{
  this->structure = structure_in;
  this->term_function = term_in;
  this->variables = variables_in;
  this->name = name_in;
  this->check = check_in;
  this->current_plate_data = data_in;
  this->structure_built = false;
  this->structure_version = 0;
}

PlateLikelihoodEstimator::~PlateLikelihoodEstimator()
{
}

PlateLikelihoodEstimator::PlateLikelihoodEstimator(const PlateLikelihoodEstimator &another)
  :ExactLikelihoodEstimator(another)
{
  this->make_copy(another);
}

void PlateLikelihoodEstimator::operator=(const PlateLikelihoodEstimator &another)
{
  if (this == &another)
    return;
  ExactLikelihoodEstimator::operator=(another);
  this->make_copy(another);
}

LikelihoodEstimator* PlateLikelihoodEstimator::duplicate() const
{
  return new PlateLikelihoodEstimator(*this);
}

void PlateLikelihoodEstimator::make_copy(const PlateLikelihoodEstimator &another)
{
  this->structure = another.structure;
  this->term_function = another.term_function;
  this->gradient_function = another.gradient_function;
  this->gradients_function = another.gradients_function;
  this->value_gradients_function = another.value_gradients_function;
  this->variables = another.variables;
  this->name = another.name;
  this->check = another.check;
  this->current_plate_data = another.current_plate_data;
  this->structure_built = another.structure_built;
  this->structure_version = another.structure_version;
  this->terms = another.terms;
  this->whole_readers = another.whole_readers;
  this->element_readers = another.element_readers;
  this->variables_read = another.variables_read;
  this->read_index = another.read_index;
}

LikelihoodEstimatorOutput* PlateLikelihoodEstimator::initialise()
{
  if (!this->structure_built)
    this->build_structure();
  return new PlateLikelihoodEstimatorOutput(this);
}

LikelihoodEstimatorOutput* PlateLikelihoodEstimator::initialise(const Parameters &parameters)
{
  return this->initialise();
}

void PlateLikelihoodEstimator::specific_change_data(Data* new_data)
{
  this->current_plate_data = new_data;
  this->build_structure();
}

const Data &PlateLikelihoodEstimator::plate_data() const
{
  static const Data no_data;
  return (this->current_plate_data != NULL) ? *this->current_plate_data : no_data;
}

double PlateLikelihoodEstimator::term(const Parameters &parameters, size_t t) const
{
  return this->term_function(parameters, this->plate_data(), this->terms, t);
}

void PlateLikelihoodEstimator::build_structure()
{
  Data structure_data = this->structure(this->plate_data());
  this->terms = structure_data["terms"];
  const arma::mat &reads = structure_data["reads"];
  size_t n_terms = this->terms.n_rows;
  size_t n_variables = this->variables.size();
  if (n_terms > size_t(UINT32_MAX))
    Rcpp::stop("The factor " + this->name + " has too many terms.");

  this->whole_readers.assign(n_variables, std::vector<uint32_t>());
  this->element_readers.assign(n_variables, std::unordered_map<uint64_t, std::vector<uint32_t>>());
  std::vector<bool> read(n_variables, false);
  for (arma::uword r = 0; r < reads.n_rows; ++r)
  {
    double t = reads(r, 0), v = reads(r, 1), row = reads(r, 2), col = reads(r, 3);
    if ((t < 0) || (size_t(t) >= n_terms) || (v < 0) || (size_t(v) >= n_variables))
      Rcpp::stop("The factor " + this->name + ": internal error in its list of dependencies.");
    size_t variable = size_t(v);
    uint32_t term_number = uint32_t(t);
    read[variable] = true;
    if ((row < 0) || (col < 0))
    {
      std::vector<uint32_t> &readers = this->whole_readers[variable];
      if (readers.empty() || readers.back() != term_number)
        readers.push_back(term_number);
    }
    else
    {
      uint64_t key = (uint64_t(row) << 32) | uint64_t(col);
      std::vector<uint32_t> &readers = this->element_readers[variable][key];
      if (readers.empty() || readers.back() != term_number)
        readers.push_back(term_number);
    }
  }
  this->variables_read.clear();
  this->read_index.clear();
  for (size_t v = 0; v < n_variables; ++v)
    if (read[v])
    {
      this->variables_read.push_back(v);
      this->read_index[this->variables[v]] = v;
    }
  this->structure_built = true;
  this->structure_version = this->structure_version + 1;
}

}
