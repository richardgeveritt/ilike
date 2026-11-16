#ifndef CHUNK_PARALLEL_H
#define CHUNK_PARALLEL_H

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <exception>
#include <functional>
#include <string>
#include <vector>
#include <memory>
#include "distributions.h"

namespace ilike
{

// Chunked generations (dev/plan.md 2.6a, G4).
//
// chunk_streams: the random-number generator of each chunk k = 0, ..., K-1:
// first (seeded as the sequential code seeds its generator) jumped k times
// (one dqrng long jump each, 2^192 steps), prepared serially, so that each
// chunk's stream depends only on the seed and k, and no two overlap.
inline std::vector<RandomNumberGenerator> chunk_streams(const RandomNumberGenerator &first,
                                                        size_t number_of_chunks)
{
  std::vector<RandomNumberGenerator> streams;
  streams.reserve(number_of_chunks);
  if (number_of_chunks == 0)
    return streams;
  streams.push_back(first);
  for (size_t k = 1; k < number_of_chunks; ++k)
  {
    std::unique_ptr<dqrng::random_64bit_generator> next = streams.back().clone(1);
    streams.push_back(*dynamic_cast<RandomNumberGenerator*>(next.get()));
  }
  return streams;
}

// Run group_body(g, first_chunk, end_chunk) for the chunks [first, end),
// split into at most number_of_groups contiguous groups, one per thread
// (each group owns its workspace and runs its chunks in order). In parallel
// (RcppParallel::parallelFor over the groups) when parallel is true and
// there is more than one group; otherwise one group, on this thread. Since
// every chunk has its own stream and writes only its own range of
// preallocated storage, the results do not depend on the grouping. An
// exception in a group is caught there and the message of the first failing
// group (in order) is raised on this thread after all groups have finished.
// The bodies must not call the R API.
struct ChunkGroupWorker : public RcppParallel::Worker
{
  const std::function<void(size_t, size_t, size_t)>* body;
  size_t first;
  size_t chunks_per_group;
  size_t end;
  std::vector<std::string>* messages;

  void operator()(std::size_t begin_group, std::size_t end_group)
  {
    for (std::size_t g = begin_group; g < end_group; ++g)
    {
      size_t a = this->first + g*this->chunks_per_group;
      size_t b = std::min(this->end, a + this->chunks_per_group);
      if (a >= b)
        continue;
      try
      {
        (*this->body)(g, a, b);
      }
      catch (std::exception &e)
      {
        (*this->messages)[g] = e.what();
        if ((*this->messages)[g].empty())
          (*this->messages)[g] = "unknown error";
      }
      catch (...)
      {
        (*this->messages)[g] = "unknown error";
      }
    }
  }
};

inline void run_chunk_groups(size_t first,
                             size_t end,
                             bool parallel,
                             size_t number_of_groups,
                             const std::function<void(size_t, size_t, size_t)> &group_body)
{
  if (end <= first)
    return;
  size_t number_of_chunks = end - first;
  if (!parallel || (number_of_groups <= 1) || (number_of_chunks <= 1))
  {
    group_body(0, first, end);
    return;
  }
  size_t groups = std::min(number_of_groups, number_of_chunks);
  size_t chunks_per_group = (number_of_chunks + groups - 1)/groups;
  groups = (number_of_chunks + chunks_per_group - 1)/chunks_per_group;
  std::vector<std::string> messages(groups);
  ChunkGroupWorker worker;
  worker.body = &group_body;
  worker.first = first;
  worker.chunks_per_group = chunks_per_group;
  worker.end = end;
  worker.messages = &messages;
  RcppParallel::parallelFor(0, groups, worker, 1);
  for (size_t g = 0; g < groups; ++g)
  {
    if (!messages[g].empty())
      Rcpp::stop(messages[g]);
  }
}

}

#endif
