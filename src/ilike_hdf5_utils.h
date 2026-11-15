#ifndef ILIKE_HDF5_UTILS_H
#define ILIKE_HDF5_UTILS_H

// Suppress HighFive's optional Boost/Eigen integrations
#ifndef HIGHFIVE_USE_BOOST
#define HIGHFIVE_USE_BOOST 0
#endif
#ifndef HIGHFIVE_USE_EIGEN
#define HIGHFIVE_USE_EIGEN 0
#endif

// Windows.h (included via RcppParallel -> tinythread.h -> windows.h) defines
// FILE_CREATE as a numeric constant, which breaks HighFive's PropertyType enum.
// Undef it and other conflicting Win32 macros before including HighFive.
#ifdef FILE_CREATE
#undef FILE_CREATE
#endif
#ifdef FILE_ACCESS
#undef FILE_ACCESS
#endif
#ifdef DATASET_CREATE
#undef DATASET_CREATE
#endif
#ifdef DATASET_ACCESS
#undef DATASET_ACCESS
#endif
#ifdef DATASET_XFER
#undef DATASET_XFER
#endif

#include <highfive/highfive.hpp>
#include <RcppArmadillo.h>
#include <string>
#include <vector>

namespace ilike
{

// Open or create an HDF5 file. Returns a shared_ptr.
inline std::shared_ptr<HighFive::File> h5_open_or_create(const std::string &path)
{
  // Build OpenOrCreate (ReadWrite|Create) from enum class values directly to
  // avoid ODR issues with the constexpr static member on macOS flat-namespace linking.
  const auto mode = HighFive::File::AccessMode::ReadWrite |
                    HighFive::File::AccessMode::Create;
  return std::make_shared<HighFive::File>(path, mode);
}

// Ensure a group exists (create if not present).
inline HighFive::Group h5_ensure_group(HighFive::File &file, const std::string &path)
{
  if (file.exist(path))
    return file.getGroup(path);
  return file.createGroup(path);
}

// Append a scalar double to a 1-D extendable dataset.
inline void h5_append_double(HighFive::File &file,
                              const std::string &dset_path,
                              double val)
{
  if (file.exist(dset_path))
  {
    auto ds = file.getDataSet(dset_path);
    auto dims = ds.getDimensions();
    size_t old_n = dims[0];
    ds.resize({old_n + 1});
    ds.select({old_n}, {1}).write(std::vector<double>{val});
  }
  else
  {
    HighFive::DataSetCreateProps props;
    props.add(HighFive::Chunking(std::vector<hsize_t>{64}));
    HighFive::DataSpace space = HighFive::DataSpace({0}, {HighFive::DataSpace::UNLIMITED});
    auto ds = file.createDataSet<double>(dset_path, space, props);
    ds.resize({1});
    ds.select({0}, {1}).write(std::vector<double>{val});
  }
}

// Append a vector of doubles to a 1-D extendable dataset in a single operation.
inline void h5_append_doubles(HighFive::File &file,
                               const std::string &dset_path,
                               const std::vector<double> &vals)
{
  if (vals.empty()) return;
  size_t n = vals.size();
  if (file.exist(dset_path))
  {
    auto ds = file.getDataSet(dset_path);
    auto dims = ds.getDimensions();
    size_t old_n = dims[0];
    ds.resize({old_n + n});
    ds.select({old_n}, {n}).write(vals);
  }
  else
  {
    HighFive::DataSetCreateProps props;
    props.add(HighFive::Chunking(std::vector<hsize_t>{std::max(n, size_t{64})}));
    HighFive::DataSpace space = HighFive::DataSpace({0}, {HighFive::DataSpace::UNLIMITED});
    auto ds = file.createDataSet<double>(dset_path, space, props);
    ds.resize({n});
    ds.select({0}, {n}).write(vals);
  }
}

// Append a scalar int to a 1-D extendable dataset.
inline void h5_append_int(HighFive::File &file,
                           const std::string &dset_path,
                           int val)
{
  if (file.exist(dset_path))
  {
    auto ds = file.getDataSet(dset_path);
    auto dims = ds.getDimensions();
    size_t old_n = dims[0];
    ds.resize({old_n + 1});
    ds.select({old_n}, {1}).write(std::vector<int>{val});
  }
  else
  {
    HighFive::DataSetCreateProps props;
    props.add(HighFive::Chunking(std::vector<hsize_t>{64}));
    HighFive::DataSpace space = HighFive::DataSpace({0}, {HighFive::DataSpace::UNLIMITED});
    auto ds = file.createDataSet<int>(dset_path, space, props);
    ds.resize({1});
    ds.select({0}, {1}).write(std::vector<int>{val});
  }
}

// Append a row to a 2-D extendable dataset (rows x cols).
// row must have the same size each call (cols is fixed by first write).
inline void h5_append_row(HighFive::File &file,
                           const std::string &dset_path,
                           const std::vector<double> &row)
{
  size_t ncols = row.size();
  if (file.exist(dset_path))
  {
    auto ds = file.getDataSet(dset_path);
    auto dims = ds.getDimensions();
    size_t old_rows = dims[0];
    ds.resize({old_rows + 1, ncols});
    ds.select({old_rows, 0}, {1, ncols}).write(std::vector<std::vector<double>>{row});
  }
  else
  {
    // Rows per chunk: HDF5 allocates whole chunks, so a chunk of 64 rows of
    // a long row (e.g. output_lengths: one row of N particles per SMC
    // generation) wasted 63/64 of the space when only one row was written
    // (IS: 512 bytes per point for N = 10^4). Keep chunks to about 8192
    // values (64 KB) while still grouping short rows.
    hsize_t rows_per_chunk = std::max<hsize_t>(1, std::min<hsize_t>(64, 8192/std::max<size_t>(ncols, 1)));
    HighFive::DataSetCreateProps props;
    props.add(HighFive::Chunking(std::vector<hsize_t>{rows_per_chunk, ncols}));
    HighFive::DataSpace space = HighFive::DataSpace({0, ncols},
                                                    {HighFive::DataSpace::UNLIMITED, ncols});
    auto ds = file.createDataSet<double>(dset_path, space, props);
    ds.resize({1, ncols});
    ds.select({0, 0}, {1, ncols}).write(std::vector<std::vector<double>>{row});
  }
}

// Append rows (given row by row in one buffer, ncols values each) to a 2-D
// extendable dataset in one operation (h5_append_row extends it once per row).
inline void h5_append_rows(HighFive::File &file,
                           const std::string &dset_path,
                           const std::vector<double> &rows_buffer,
                           size_t ncols)
{
  if ((ncols == 0) || rows_buffer.empty()) return;
  size_t nrows = rows_buffer.size()/ncols;
  if (file.exist(dset_path))
  {
    auto ds = file.getDataSet(dset_path);
    size_t old_rows = ds.getDimensions()[0];
    ds.resize({old_rows + nrows, ncols});
    ds.select({old_rows, 0}, {nrows, ncols}).write_raw(rows_buffer.data());
  }
  else
  {
    hsize_t rows_per_chunk = std::max<hsize_t>(1, std::min<hsize_t>(std::max<size_t>(nrows, 64), 8192/ncols));
    HighFive::DataSetCreateProps props;
    props.add(HighFive::Chunking(std::vector<hsize_t>{rows_per_chunk, ncols}));
    HighFive::DataSpace space = HighFive::DataSpace({nrows, ncols},
                                                    {HighFive::DataSpace::UNLIMITED, ncols});
    auto ds = file.createDataSet<double>(dset_path, space, props);
    ds.write_raw(rows_buffer.data());
  }
}

// As h5_append_rows, from a raw buffer of nrows rows of ncols values (row by
// row): e.g. the memory of a (ncols x nrows) Armadillo matrix.
inline void h5_append_rows(HighFive::File &file,
                           const std::string &dset_path,
                           const double *rows_buffer,
                           size_t nrows,
                           size_t ncols)
{
  if ((ncols == 0) || (nrows == 0)) return;
  if (file.exist(dset_path))
  {
    auto ds = file.getDataSet(dset_path);
    size_t old_rows = ds.getDimensions()[0];
    ds.resize({old_rows + nrows, ncols});
    ds.select({old_rows, 0}, {nrows, ncols}).write_raw(rows_buffer);
  }
  else
  {
    hsize_t rows_per_chunk = std::max<hsize_t>(1, std::min<hsize_t>(std::max<size_t>(nrows, 64), 8192/ncols));
    HighFive::DataSetCreateProps props;
    props.add(HighFive::Chunking(std::vector<hsize_t>{rows_per_chunk, ncols}));
    HighFive::DataSpace space = HighFive::DataSpace({nrows, ncols},
                                                    {HighFive::DataSpace::UNLIMITED, ncols});
    auto ds = file.createDataSet<double>(dset_path, space, props);
    ds.write_raw(rows_buffer);
  }
}

// Write a vector<string> as a fixed string attribute on a group or dataset.
inline void h5_set_str_attr(HighFive::Group &grp,
                             const std::string &attr_name,
                             const std::vector<std::string> &vals)
{
  if (grp.hasAttribute(attr_name))
    grp.deleteAttribute(attr_name);
  grp.createAttribute(attr_name, vals);
}

// Write a vector<size_t> as a fixed attribute on a group.
inline void h5_set_sizet_attr(HighFive::Group &grp,
                               const std::string &attr_name,
                               const std::vector<size_t> &vals)
{
  if (grp.hasAttribute(attr_name))
    grp.deleteAttribute(attr_name);
  grp.createAttribute(attr_name, vals);
}

// Write an arma::mat as a 2-D dataset (rows = observations, cols = variables).
// Stored row-major: element (i,j) is row i, col j.
inline void h5_write_mat(HighFive::Group &grp,
                          const std::string &name,
                          const arma::mat &m)
{
  size_t nrows = m.n_rows;
  size_t ncols = m.n_cols;
  // Armadillo is column-major; copy to row-major std::vector
  std::vector<std::vector<double>> data(nrows, std::vector<double>(ncols));
  for (size_t r = 0; r < nrows; ++r)
    for (size_t c = 0; c < ncols; ++c)
      data[r][c] = m(r, c);

  if (grp.exist(name))
    grp.unlink(name);
  grp.createDataSet(name, data);
}

// Write a (columns x points) matrix as a (points x columns) dataset: the
// column-major memory of the argument is exactly the row-major layout of the
// dataset, so it is written in one call without any copy (h5_write_mat copies
// every row into its own std::vector first).
inline void h5_write_mat_transposed(HighFive::Group &grp,
                                    const std::string &name,
                                    const arma::mat &transposed)
{
  if (grp.exist(name))
    grp.unlink(name);
  HighFive::DataSpace space(std::vector<size_t>{transposed.n_cols, transposed.n_rows});
  auto ds = grp.createDataSet<double>(name, space);
  ds.write_raw(transposed.memptr());
}

// Write an arma::colvec as a 1-D dataset.
inline void h5_write_vec(HighFive::Group &grp,
                          const std::string &name,
                          const arma::colvec &v)
{
  std::vector<double> data(v.begin(), v.end());
  if (grp.exist(name))
    grp.unlink(name);
  grp.createDataSet(name, data);
}

// Write a vector<size_t> as a 1-D dataset.
inline void h5_write_intvec(HighFive::Group &grp,
                              const std::string &name,
                              const std::vector<size_t> &v)
{
  if (grp.exist(name))
    grp.unlink(name);
  grp.createDataSet(name, v);
}

// Write a scalar double as a dataset.
inline void h5_write_scalar_double(HighFive::Group &grp,
                                    const std::string &name,
                                    double val)
{
  if (grp.exist(name))
    grp.unlink(name);
  grp.createDataSet(name, val);
}

// Write a scalar int as a dataset.
inline void h5_write_scalar_int(HighFive::Group &grp,
                                 const std::string &name,
                                 int val)
{
  if (grp.exist(name))
    grp.unlink(name);
  grp.createDataSet(name, val);
}

// Write a string as a dataset.
inline void h5_write_string(HighFive::Group &grp,
                              const std::string &name,
                              const std::string &val)
{
  if (grp.exist(name))
    grp.unlink(name);
  grp.createDataSet(name, val);
}

} // namespace ilike

#endif // ILIKE_HDF5_UTILS_H
