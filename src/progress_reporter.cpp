#include <RcppArmadillo.h>
#include <cli/progress.h>
#include "progress_reporter.h"
#include "chunk_parallel.h"

namespace ilike
{

ProgressReporter& ProgressReporter::instance()
{
  static ProgressReporter reporter;
  return reporter;
}

static bool on_main_thread_outside_parallel()
{
  return (std::this_thread::get_id() == ilike_main_thread()) && (ilike_parallel_depth() == 0);
}

void ProgressReporter::begin(double total, const std::string &name)
{
  if (this->active || !on_main_thread_outside_parallel())
    return;
  SEXP config = PROTECT(Rf_allocVector(VECSXP, 1));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 1));
  SET_VECTOR_ELT(config, 0, Rf_mkString(name.c_str()));
  SET_STRING_ELT(names, 0, Rf_mkChar("name"));
  Rf_setAttrib(config, R_NamesSymbol, names);
  SEXP the_bar = cli_progress_bar((total < 0.0) ? NA_REAL : total, config);
  R_PreserveObject(the_bar);
  UNPROTECT(2);
  this->bar = (void*)the_bar;
  this->active = true;
}

void ProgressReporter::update(double position, const std::string &status)
{
  if (!this->active || !CLI_SHOULD_TICK || !on_main_thread_outside_parallel())
    return;
  SEXP the_bar = (SEXP)this->bar;
  if (!status.empty())
    cli_progress_set_status(the_bar, status.c_str());
  cli_progress_set(the_bar, position);
}

void ProgressReporter::end()
{
  if (!this->active)
    return;
  SEXP the_bar = (SEXP)this->bar;
  if (on_main_thread_outside_parallel())
    cli_progress_done(the_bar);
  R_ReleaseObject(the_bar);
  this->bar = nullptr;
  this->active = false;
}

}
