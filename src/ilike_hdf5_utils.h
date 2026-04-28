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
    HighFive::DataSetCreateProps props;
    props.add(HighFive::Chunking(std::vector<hsize_t>{64, ncols}));
    HighFive::DataSpace space = HighFive::DataSpace({0, ncols},
                                                    {HighFive::DataSpace::UNLIMITED, ncols});
    auto ds = file.createDataSet<double>(dset_path, space, props);
    ds.resize({1, ncols});
    ds.select({0, 0}, {1, ncols}).write(std::vector<std::vector<double>>{row});
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
