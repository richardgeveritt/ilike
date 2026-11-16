#ifndef OPTIONAL_DELETE_H
#define OPTIONAL_DELETE_H

namespace ilike
{

// A deleter for std::unique_ptr when the pointer is owned only sometimes:
// std::unique_ptr<T, OptionalDelete> deletes its object only if it was given
// OptionalDelete(true). Used by the MCMC outputs, which own their MCMC move
// when the move adapts (each chain adapts its own copy) and otherwise share
// the algorithm's move without owning it. Whether the pointer is owned is
// then part of the pointer (get_deleter().owns), not a separate flag.
struct OptionalDelete
{
  bool owns = true;

  OptionalDelete() = default;
  explicit OptionalDelete(bool owns_in) : owns(owns_in) {}

  template <class T>
  void operator()(T *pointer) const
  {
    if (this->owns)
      delete pointer;
  }
};

}

#endif
