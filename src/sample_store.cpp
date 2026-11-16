#include "sample_store.h"

namespace ilike
{

SampleStore::SampleStore()
: layout_set(false), total_dimension(0), number_of_extra_columns(0), number_of_points(0)
{
}

SampleStore::~SampleStore()
{
}

SampleStore::SampleStore(const SampleStore &another) = default;
SampleStore& SampleStore::operator=(const SampleStore &another) = default;
SampleStore::SampleStore(SampleStore &&another) = default;
SampleStore& SampleStore::operator=(SampleStore &&another) = default;

void SampleStore::set_number_of_extra_columns(size_t number_of_extra_columns_in)
{
  if (this->number_of_points > 0)
    Rcpp::stop("SampleStore::set_number_of_extra_columns: the store already holds points.");
  this->number_of_extra_columns = number_of_extra_columns_in;
}

size_t SampleStore::get_number_of_extra_columns() const
{
  return this->number_of_extra_columns;
}

void SampleStore::setup_layout(const Parameters &parameters)
{
  this->vector_names = parameters.get_vector_variables();
  this->vector_offsets.clear();
  this->vector_n_rows.clear();
  this->vector_n_cols.clear();
  this->total_dimension = 0;
  for (auto i = this->vector_names.begin(); i != this->vector_names.end(); ++i)
  {
    const arma::mat &m = parameters[*i];
    this->vector_offsets.push_back(this->total_dimension);
    this->vector_n_rows.push_back(m.n_rows);
    this->vector_n_cols.push_back(m.n_cols);
    this->total_dimension += m.n_elem;
  }
  this->any_names = parameters.get_any_variables();
  this->any_values.assign(this->any_names.size(), std::vector<boost::any>());
  this->layout_set = true;
}

bool SampleStore::compatible(const Parameters &parameters) const
{
  if (!this->layout_set)
    return true;
  std::vector<std::string> names = parameters.get_vector_variables();
  if (names.size() != this->vector_names.size())
    return false;
  for (size_t k = 0; k < this->vector_names.size(); ++k)
  {
    if (std::find(names.begin(), names.end(), this->vector_names[k]) == names.end())
      return false;
    const arma::mat &m = parameters[this->vector_names[k]];
    if ((m.n_rows != this->vector_n_rows[k]) || (m.n_cols != this->vector_n_cols[k]))
      return false;
  }
  std::vector<std::string> any_names_in = parameters.get_any_variables();
  if (any_names_in.size() != this->any_names.size())
    return false;
  for (size_t k = 0; k < this->any_names.size(); ++k)
  {
    if (std::find(any_names_in.begin(), any_names_in.end(), this->any_names[k]) == any_names_in.end())
      return false;
  }
  return true;
}

void SampleStore::grow(size_t minimum_capacity)
{
  size_t capacity = this->values.n_cols;
  if (minimum_capacity <= capacity)
    return;
  size_t new_capacity = std::max(minimum_capacity, std::max<size_t>(2*capacity, 16));
  this->values.resize(this->total_dimension, new_capacity);
  this->extras.resize(this->number_of_extra_columns, new_capacity);
}

const arma::mat& SampleStore::find_vector_variable(const Parameters &parameters,
                                                   size_t k)
{
  if (this->vector_positions.size() != this->vector_names.size())
    this->vector_positions.assign(this->vector_names.size(), 0);
  const std::string &name = this->vector_names[k];
  size_t position = this->vector_positions[k];
  if (position < parameters.vector_size())
  {
    auto i = parameters.vector_begin() + position;
    if (i->first == name)
      return *i->second.first;
  }
  size_t j = 0;
  for (auto i = parameters.vector_begin(); i != parameters.vector_end(); ++i, ++j)
  {
    if (i->first == name)
    {
      this->vector_positions[k] = j;
      return *i->second.first;
    }
  }
  // as before: a missing variable stops
  return parameters[name];
}

void SampleStore::append(const Parameters &parameters,
                         const double *extra_values)
{
  // The layout is checked variable by variable below (a missing variable
  // stops in Parameters::operator[], a change of shape stops here); the full
  // comparison of variable names (compatible) is not repeated for every
  // point, since building the list of names costs more than storing the point.
  if (!this->layout_set)
    this->setup_layout(parameters);

  if (this->pending_reserve > 0)
  {
    this->grow(this->pending_reserve);
    this->pending_reserve = 0;
  }
  this->grow(this->number_of_points + 1);

  double *column = this->values.colptr(this->number_of_points);
  for (size_t k = 0; k < this->vector_names.size(); ++k)
  {
    // element order of as_row(): row by row
    const arma::mat &m = this->find_vector_variable(parameters, k);
    if ((m.n_rows != this->vector_n_rows[k]) || (m.n_cols != this->vector_n_cols[k]))
      Rcpp::stop("SampleStore::append: variable '" + this->vector_names[k] + "' changed size.");
    size_t offset = this->vector_offsets[k];
    for (arma::uword r = 0; r < m.n_rows; ++r)
      for (arma::uword c = 0; c < m.n_cols; ++c)
        column[offset++] = m(r, c);
  }
  for (size_t k = 0; k < this->any_names.size(); ++k)
  {
    this->any_values[k].push_back(parameters(this->any_names[k]));
  }
  if (this->number_of_extra_columns > 0)
  {
    double *extra_column = this->extras.colptr(this->number_of_points);
    for (size_t j = 0; j < this->number_of_extra_columns; ++j)
      extra_column[j] = (extra_values == NULL) ? arma::datum::nan : extra_values[j];
  }
  this->number_of_points = this->number_of_points + 1;
}

const arma::mat& SampleStore::lookup_vector_variable(const Parameters &parameters,
                                                     size_t k) const
{
  const std::string &name = this->vector_names[k];
  if (k < parameters.vector_size())
  {
    auto i = parameters.vector_begin() + k;
    if (i->first == name)
      return *i->second.first;
  }
  return parameters[name];
}

void SampleStore::prepare(const Parameters &parameters, size_t number_of_points_in)
{
  if (!this->layout_set)
    this->setup_layout(parameters);
  this->pending_reserve = 0;
  this->grow(number_of_points_in);
  for (auto i = this->any_values.begin(); i != this->any_values.end(); ++i)
    i->resize(number_of_points_in);
  this->number_of_points = number_of_points_in;
}

void SampleStore::set_point(size_t i,
                            const Parameters &parameters,
                            const double *extra_values)
{
  if (i >= this->number_of_points)
    Rcpp::stop("SampleStore::set_point: index out of range.");
  double *column = this->values.colptr(i);
  for (size_t k = 0; k < this->vector_names.size(); ++k)
  {
    // element order of as_row(): row by row
    const arma::mat &m = this->lookup_vector_variable(parameters, k);
    if ((m.n_rows != this->vector_n_rows[k]) || (m.n_cols != this->vector_n_cols[k]))
      Rcpp::stop("SampleStore::set_point: variable '" + this->vector_names[k] + "' changed size.");
    size_t offset = this->vector_offsets[k];
    for (arma::uword r = 0; r < m.n_rows; ++r)
      for (arma::uword c = 0; c < m.n_cols; ++c)
        column[offset++] = m(r, c);
  }
  for (size_t k = 0; k < this->any_names.size(); ++k)
  {
    this->any_values[k][i] = parameters(this->any_names[k]);
  }
  if (this->number_of_extra_columns > 0)
  {
    double *extra_column = this->extras.colptr(i);
    for (size_t j = 0; j < this->number_of_extra_columns; ++j)
      extra_column[j] = (extra_values == NULL) ? arma::datum::nan : extra_values[j];
  }
}

size_t SampleStore::size() const
{
  return this->number_of_points;
}

bool SampleStore::empty() const
{
  return this->number_of_points == 0;
}

void SampleStore::reserve(size_t number_of_points_in)
{
  if (this->layout_set)
    this->grow(number_of_points_in);
}

void SampleStore::clear()
{
  this->number_of_points = 0;
  for (auto i = this->any_values.begin(); i != this->any_values.end(); ++i)
    i->clear();
}

arma::mat SampleStore::matrix_of(const std::vector<std::string> &variables) const
{
  std::vector<size_t> which;
  size_t number_of_columns = 0;
  for (auto v = variables.begin(); v != variables.end(); ++v)
  {
    auto found = std::find(this->vector_names.begin(), this->vector_names.end(), *v);
    if (found == this->vector_names.end())
      Rcpp::stop("SampleStore::matrix_of: variable '" + *v + "' not found in the stored points.");
    size_t k = size_t(found - this->vector_names.begin());
    which.push_back(k);
    number_of_columns += this->vector_n_rows[k]*this->vector_n_cols[k];
  }
  arma::mat output(this->number_of_points, number_of_columns);
  size_t output_column = 0;
  for (auto k = which.begin(); k != which.end(); ++k)
  {
    size_t n = this->vector_n_rows[*k]*this->vector_n_cols[*k];
    if (this->number_of_points > 0)
      output.cols(output_column, output_column + n - 1) =
        this->values.submat(this->vector_offsets[*k], 0,
                            this->vector_offsets[*k] + n - 1, this->number_of_points - 1).t();
    output_column += n;
  }
  return output;
}

void SampleStore::fill_transposed(const std::vector<std::string> &variables,
                                  arma::mat &destination,
                                  size_t first_column) const
{
  if (this->number_of_points == 0)
    return;
  size_t destination_row = 0;
  for (auto v = variables.begin(); v != variables.end(); ++v)
  {
    auto found = std::find(this->vector_names.begin(), this->vector_names.end(), *v);
    if (found == this->vector_names.end())
      Rcpp::stop("SampleStore::fill_transposed: variable '" + *v + "' not found in the stored points.");
    size_t k = size_t(found - this->vector_names.begin());
    size_t n = this->vector_n_rows[k]*this->vector_n_cols[k];
    destination.submat(destination_row, first_column,
                       destination_row + n - 1, first_column + this->number_of_points - 1) =
      this->values.submat(this->vector_offsets[k], 0,
                          this->vector_offsets[k] + n - 1, this->number_of_points - 1);
    destination_row += n;
  }
}

void SampleStore::reserve_after_first(size_t number_of_points_in)
{
  if (this->layout_set)
    this->grow(number_of_points_in);
  else
    this->pending_reserve = number_of_points_in;
}

Parameters SampleStore::parameters_at(size_t i) const
{
  if (i >= this->number_of_points)
    Rcpp::stop("SampleStore::parameters_at: index out of range.");
  Parameters output;
  const double *column = this->values.colptr(i);
  for (size_t k = 0; k < this->vector_names.size(); ++k)
  {
    arma::mat m(this->vector_n_rows[k], this->vector_n_cols[k]);
    size_t offset = this->vector_offsets[k];
    for (arma::uword r = 0; r < m.n_rows; ++r)
      for (arma::uword c = 0; c < m.n_cols; ++c)
        m(r, c) = column[offset++];
    output[this->vector_names[k]] = m;
  }
  for (size_t k = 0; k < this->any_names.size(); ++k)
  {
    output(this->any_names[k]) = this->any_values[k][i];
  }
  return output;
}

void SampleStore::assign_into(size_t i, Parameters &target) const
{
  if (i >= this->number_of_points)
    Rcpp::stop("SampleStore::assign_into: index out of range.");
  const double *column = this->values.colptr(i);
  for (size_t k = 0; k < this->vector_names.size(); ++k)
  {
    const std::string &name = this->vector_names[k];
    arma::uword n_rows = this->vector_n_rows[k];
    arma::uword n_cols = this->vector_n_cols[k];
    arma::mat* destination = NULL;
    bool found = false;
    for (auto v = target.vector_begin(); v != target.vector_end(); ++v)
    {
      if (v->first != name)
        continue;
      found = true;
      if (v->second.second)
        break; // fixed: not ours to write
      if ((v->second.first.use_count() != 1) ||
          (v->second.first->n_rows != n_rows) || (v->second.first->n_cols != n_cols))
        v->second.first = std::make_shared<arma::mat>(n_rows, n_cols);
      destination = v->second.first.get();
      break;
    }
    if (!found)
    {
      target[name] = arma::mat(n_rows, n_cols);
      destination = &target[name];
    }
    if (destination == NULL)
      continue;
    // element order of as_row(): row by row
    size_t offset = this->vector_offsets[k];
    for (arma::uword r = 0; r < n_rows; ++r)
      for (arma::uword c = 0; c < n_cols; ++c)
        (*destination)(r, c) = column[offset++];
  }
  for (size_t k = 0; k < this->any_names.size(); ++k)
  {
    target(this->any_names[k]) = this->any_values[k][i];
  }
}

void SampleStore::set_extras(size_t i, const double *extra_values)
{
  if (i >= this->number_of_points)
    Rcpp::stop("SampleStore::set_extras: index out of range.");
  double *extra_column = this->extras.colptr(i);
  for (size_t j = 0; j < this->number_of_extra_columns; ++j)
    extra_column[j] = extra_values[j];
}

double SampleStore::extra_at(size_t i, size_t j) const
{
  if ((i >= this->number_of_points) || (j >= this->number_of_extra_columns))
    Rcpp::stop("SampleStore::extra_at: index out of range.");
  return this->extras(j, i);
}

const double* SampleStore::extras_of(size_t i) const
{
  if ((i >= this->number_of_points) || (this->number_of_extra_columns == 0))
    Rcpp::stop("SampleStore::extras_of: index out of range.");
  return this->extras.colptr(i);
}

size_t SampleStore::memory_bytes() const
{
  return sizeof(double)*(this->values.n_elem + this->extras.n_elem);
}

}
