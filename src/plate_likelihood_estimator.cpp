#include <unordered_set>
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
  this->tempered = another.tempered;
  this->tempering_variable = another.tempering_variable;
  this->tempering_anneal = another.tempering_anneal;
  this->tempering_random = another.tempering_random;
  this->tempering_seed = another.tempering_seed;
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

void PlateLikelihoodEstimator::set_data_tempering(const std::string &variable_in, bool anneal_in, bool random_order_in, size_t seed_in)
{
  this->tempered = true;
  this->tempering_variable = variable_in;
  this->tempering_anneal = anneal_in;
  this->tempering_random = random_order_in;
  this->tempering_seed = seed_in;
  if (this->structure_built)
    this->build_structure();
}

void PlateLikelihoodEstimator::tempering_point(const Parameters &parameters, size_t &k, double &r) const
{
  size_t n = this->number_of_terms();
  double phi = parameters[this->tempering_variable](0);
  if (!(phi > 0.0))
    phi = 0.0;
  if (phi > 1.0)
    phi = 1.0;
  double x = double(n) * phi;
  r = 0.0;
  if (this->tempering_anneal)
  {
    k = size_t(std::floor(x));
    if (k >= n)
      k = n;
    else
      r = x - double(k);
  }
  else
  {
    k = size_t(std::ceil(x));
    if (k > n)
      k = n;
  }
}

// a permutation of 0, ..., n - 1 from a seed (Fisher-Yates with splitmix64,
// so that the order is the same on every platform)
static std::vector<size_t> tempering_permutation(size_t n, size_t seed)
{
  std::vector<size_t> order(n);
  for (size_t i = 0; i < n; ++i)
    order[i] = i;
  uint64_t state = uint64_t(seed) + 0x9E3779B97F4A7C15ULL;
  auto next = [&state]() {
    uint64_t z = (state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
  };
  for (size_t i = n; i > 1; --i)
  {
    size_t j = size_t(next() % uint64_t(i));
    std::swap(order[i - 1], order[j]);
  }
  return order;
}

void PlateLikelihoodEstimator::build_structure()
{
  Data structure_data = this->structure(this->plate_data());
  this->terms = structure_data["terms"];
  arma::mat reads = structure_data["reads"];
  size_t n_terms = this->terms.n_rows;
  if (this->tempered && this->tempering_random && (n_terms > 1))
  {
    // (data tempering) term j of the plate is term order[j] of the
    // statements: the terms' rows, and the term numbers in the reads, follow
    std::vector<size_t> order = tempering_permutation(n_terms, this->tempering_seed);
    std::vector<size_t> position(n_terms);
    arma::mat permuted(this->terms.n_rows, this->terms.n_cols);
    for (size_t j = 0; j < n_terms; ++j)
    {
      permuted.row(j) = this->terms.row(order[j]);
      position[order[j]] = j;
    }
    this->terms = permuted;
    for (arma::uword r = 0; r < reads.n_rows; ++r)
      if ((reads(r, 0) >= 0) && (size_t(reads(r, 0)) < n_terms))
        reads(r, 0) = double(position[size_t(reads(r, 0))]);
    // (the reads grouped by term, as below expects)
    arma::uvec by_term = arma::stable_sort_index(reads.col(0));
    reads = reads.rows(by_term);
  }
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


// (dev/plan.md 3.22 U5b, C312) from the records of what each term reads
int PlateLikelihoodEstimator::terms_read_one_element(const std::string &variable) const
{
  if (!this->structure_built)
    return -1;
  auto found = this->read_index.find(variable);
  if (found == this->read_index.end())
    return 1;
  const size_t v = found->second;
  if ((v < this->whole_readers.size()) && !this->whole_readers[v].empty())
    return 0;
  if (v >= this->element_readers.size())
    return 1;
  std::unordered_set<uint32_t> seen;
  for (const auto &element : this->element_readers[v])
    for (uint32_t term : element.second)
      if (!seen.insert(term).second)
        return 0;
  return 1;
}

}
