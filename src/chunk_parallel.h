#ifndef CHUNK_PARALLEL_H
#define CHUNK_PARALLEL_H

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <exception>
#include <functional>
#include <string>
#include <vector>
#include <memory>
#include <thread>
#include <cfenv>
#include <atomic>
#include "distributions.h"

namespace ilike
{

// Chunked generations (dev/plan.md 2.6a, G4).
//
// User interrupts (dev/plan.md 3.11). R may only be asked from its main
// thread, and not while worker threads are running (a parallel region):
// check_interrupt() does nothing elsewhere. The main thread is recorded when
// the package is loaded (ilike_record_main_thread, an Rcpp init hook). An
// interrupt unwinds the C++ stack (Rcpp's InterruptedException).
inline std::thread::id& ilike_main_thread()
{
  static std::thread::id id;
  return id;
}

inline std::atomic<int>& ilike_parallel_depth()
{
  static std::atomic<int> depth(0);
  return depth;
}

inline void check_interrupt()
{
  if ((std::this_thread::get_id() == ilike_main_thread()) && (ilike_parallel_depth() == 0))
  {
    Rcpp::checkUserInterrupt();
  }
}

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
  // The floating-point environment of the thread that started the region
  // (rounding, and on Windows/MinGW the x87 precision, which worker threads
  // do not inherit: the maths library's exp and log then differed in the
  // last bits between threads, so parallel runs differed from sequential
  // ones in the 16th digit). Each group runs under it.
  std::fenv_t environment;
  const std::function<void(size_t, size_t, size_t)>* body;
  size_t first;
  size_t chunks_per_group;
  size_t end;
  std::vector<std::string>* messages;

  void operator()(std::size_t begin_group, std::size_t end_group)
  {
    std::fenv_t own_environment;
    std::fegetenv(&own_environment);
    std::fesetenv(&this->environment);
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
    std::fesetenv(&own_environment);
  }
};

inline void run_chunk_groups(size_t first,
                             size_t end,
                             bool parallel,
                             size_t number_of_groups,
                             const std::function<void(size_t, size_t, size_t)> &group_body,
                             const std::function<void(size_t)> &after_wave = std::function<void(size_t)>())
{
  // The chunks are processed in waves, with a check for a user interrupt
  // after each (on the main thread, between parallel regions): 16 chunks per
  // wave sequentially, 8 per group in parallel. A wave's chunks are split
  // into contiguous groups as before; the grouping does not change results.
  if (end <= first)
    return;
  bool in_parallel = parallel && (number_of_groups > 1) && (end - first > 1);
  if (!in_parallel)
  {
    const size_t wave = 16;
    for (size_t a = first; a < end; a += wave)
    {
      group_body(0, a, std::min(end, a + wave));
      if (after_wave)
        after_wave(std::min(end, a + wave));
      check_interrupt();
    }
    return;
  }
  size_t wave = 8*number_of_groups;
  for (size_t wave_first = first; wave_first < end; wave_first += wave)
  {
    size_t wave_end = std::min(end, wave_first + wave);
    size_t number_of_chunks = wave_end - wave_first;
    size_t groups = std::min(number_of_groups, number_of_chunks);
    size_t chunks_per_group = (number_of_chunks + groups - 1)/groups;
    groups = (number_of_chunks + chunks_per_group - 1)/chunks_per_group;
    std::vector<std::string> messages(groups);
    ChunkGroupWorker worker;
    worker.body = &group_body;
    worker.first = wave_first;
    worker.chunks_per_group = chunks_per_group;
    worker.end = wave_end;
    worker.messages = &messages;
    std::fegetenv(&worker.environment);
    ilike_parallel_depth()++;
    RcppParallel::parallelFor(0, groups, worker, 1);
    ilike_parallel_depth()--;
    for (size_t g = 0; g < groups; ++g)
    {
      if (!messages[g].empty())
        Rcpp::stop(messages[g]);
    }
    if (after_wave)
      after_wave(wave_end);
    check_interrupt();
  }
}

}

#endif
