#ifndef PROPOSALTARGET_H
#define PROPOSALTARGET_H

#include <RcppArmadillo.h>
#include <cstdint>
#include "parameters.h"

namespace ilike
{

// The parameters of a proposal being built in a scratch particle
// (dev/plan.md, Stage 3 item 6). The result is always that of
//   target = base; target.deep_overwrite_with_variables_in_argument(proposed);
// where `proposed` holds the variables written through matrix() / write_all().
// When the target already has base's variables (same names, same order: every
// proposal after the first), a written variable goes straight into the
// target's own matrix for it — reused when no other Parameters object shares
// it — and every variable not written shares base's matrix, so a kernel that
// writes through matrix() proposes without allocating. Otherwise (first use,
// a new variable, non-matrix variables) the written values are collected and
// the two steps above are done literally.
class ProposalTarget
{
public:

  ProposalTarget(Parameters &target_in,
                 const Parameters &base_in)
  : target(target_in), base(base_in), written(0), collect(false)
  {
    this->direct = (this->target.vector_size() == this->base.vector_size()) &&
                   (this->target.vector_size() <= 64) &&
                   (this->target.any_size() == 0) && (this->base.any_size() == 0);
    if (this->direct)
    {
      auto j = this->target.vector_begin();
      for (auto i = this->base.vector_begin(); i != this->base.vector_end(); ++i, ++j)
      {
        if (i->first != j->first)
        {
          this->direct = false;
          break;
        }
      }
    }
  }

  // The matrix to write the proposed value of `name` into (assign the whole
  // value: its previous contents are stale).
  arma::mat& matrix(const std::string &name)
  {
    if (this->direct && !this->collect)
    {
      size_t k = 0;
      auto b = this->base.vector_begin();
      for (auto j = this->target.vector_begin(); j != this->target.vector_end(); ++j, ++b, ++k)
      {
        if (j->first == name)
        {
          if (b->second.second == true)
            Rcpp::stop("Parameters::deep_overwrite_with_variables_in_argument - parameter is fixed and cannot be changed.");
          if (j->second.first.use_count() != 1)
            j->second.first = std::make_shared<arma::mat>();
          j->second.second = false;
          this->written |= (uint64_t(1) << k);
          return *j->second.first;
        }
      }
      // a variable base does not have: fall back to the literal two steps
      this->start_collecting();
    }
    return this->collected[name];
  }

  // Write every variable of `values` (the default for kernels that return a
  // new Parameters from simulate).
  void write_all(const Parameters &values)
  {
    for (auto i = values.vector_begin(); i != values.vector_end(); ++i)
      this->matrix(i->first) = *i->second.first;
    if (values.any_size() > 0)
    {
      if (this->direct)
        this->start_collecting();
      for (auto i = values.any_begin(); i != values.any_end(); ++i)
        this->collected(i->first) = *i->second.first;
    }
  }

  // Complete the target: call once, after all writes.
  void finish()
  {
    if (this->direct && !this->collect)
    {
      size_t k = 0;
      auto i = this->base.vector_begin();
      for (auto j = this->target.vector_begin(); j != this->target.vector_end(); ++j, ++i, ++k)
      {
        if (!(this->written & (uint64_t(1) << k)))
          j->second = i->second;
      }
      return;
    }
    this->target = this->base;
    this->target.deep_overwrite_with_variables_in_argument(this->collected);
  }

private:

  // move what was written in place so far into `collected`
  void start_collecting()
  {
    this->collect = true;
    size_t k = 0;
    for (auto j = this->target.vector_begin(); j != this->target.vector_end(); ++j, ++k)
    {
      if (this->written & (uint64_t(1) << k))
        this->collected[j->first] = *j->second.first;
    }
  }

  Parameters &target;
  const Parameters &base;
  bool direct;          // target has base's layout: write in place
  uint64_t written;     // positions written in place
  bool collect;         // gave up writing in place
  Parameters collected; // the written values, when not in place
};

}

#endif
