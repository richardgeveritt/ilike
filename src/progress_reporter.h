#ifndef PROGRESS_REPORTER_H
#define PROGRESS_REPORTER_H

#include <string>

namespace ilike
{

// Progress bars (dev/plan.md 3.11), through cli's C API. One reporter for
// the whole package (the run that asked for one owns it; inner algorithms
// never start one). Updates go to cli only on R's main thread, outside
// parallel regions, and only when cli's timer says an update is due, so a
// call costs a branch when nothing is shown. All cli calls are in
// progress_reporter.cpp (cli's timer flag is per compilation unit).
class ProgressReporter
{
public:
  static ProgressReporter& instance();
  // total < 0: unknown (a spinner)
  void begin(double total, const std::string &name);
  void update(double position, const std::string &status = "");
  void end();
  bool is_active() const { return this->active; }
private:
  bool active = false;
  void* bar = nullptr; // the cli bar (an R object, preserved while active)
};

// Ends the bar when it goes out of scope (also on an error or interrupt).
struct ProgressScope
{
  bool owns;
  ProgressScope(bool show, double total, const std::string &name) : owns(show)
  {
    if (owns)
      ProgressReporter::instance().begin(total, name);
  }
  ~ProgressScope()
  {
    if (owns)
      ProgressReporter::instance().end();
  }
};

}

#endif
