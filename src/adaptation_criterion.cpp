#include "adaptation_criterion.h"
#include "mcmc_adaptor.h"

namespace ilike
{

AdaptationCriterion::AdaptationCriterion()
: kind(Kind::iterations), n(0), tolerance(0.0), window(1), relative(true)
{
}

AdaptationCriterion::AdaptationCriterion(const AdaptationCriterion &another)
{
  this->make_copy(another);
}

AdaptationCriterion &AdaptationCriterion::operator=(const AdaptationCriterion &another)
{
  if (this != &another)
    this->make_copy(another);
  return *this;
}

AdaptationCriterion* AdaptationCriterion::duplicate() const
{
  return new AdaptationCriterion(*this);
}

void AdaptationCriterion::make_copy(const AdaptationCriterion &another)
{
  this->kind = another.kind;
  this->n = another.n;
  this->tolerance = another.tolerance;
  this->window = another.window;
  this->watch = another.watch;
  this->relative = another.relative;
  this->history = another.history;
  this->acceptances = another.acceptances;
  this->items.clear();
  for (const auto &item : another.items)
    this->items.push_back(std::unique_ptr<AdaptationCriterion>(item->duplicate()));
}

AdaptationCriterion* AdaptationCriterion::from_recipe(const RecipeValue &node)
{
  std::unique_ptr<AdaptationCriterion> output(new AdaptationCriterion());
  std::string criterion = node["criterion"].as_string();
  if (criterion == "iterations")
  {
    output->kind = Kind::iterations;
    output->n = node["n"].as_size_t();
  }
  else if ((criterion == "stable") || (criterion == "acceptance"))
  {
    output->kind = (criterion == "stable") ? Kind::stable : Kind::acceptance;
    output->tolerance = node["tolerance"].as_double();
    output->window = node["window"].as_size_t();
    if (output->window < 1)
      Rcpp::stop(criterion + ": the window is at least 1.");
    if (node.contains("watch"))
      output->watch = node["watch"].as_strings();
    if (node.contains("change"))
      output->relative = (node["change"].as_string() != "absolute");
  }
  else if ((criterion == "all") || (criterion == "any"))
  {
    output->kind = (criterion == "all") ? Kind::all : Kind::any;
    const RecipeValue &items = node["items"];
    for (size_t k = 0; k < items.size(); ++k)
      output->items.push_back(std::unique_ptr<AdaptationCriterion>(from_recipe(items[k])));
  }
  else
    Rcpp::stop("Unknown criterion to stop adapting: " + criterion + ".");
  return output.release();
}

bool AdaptationCriterion::met(const MCMCAdaptor &adaptor, double acceptance_probability, size_t steps)
{
  switch (this->kind)
  {
    case Kind::iterations:
      return steps >= this->n;
    case Kind::all:
    case Kind::any:
    {
      // every criterion inside sees every step (their windows stay current)
      bool all = true, any = false;
      for (auto &item : this->items)
      {
        bool one = item->met(adaptor, acceptance_probability, steps);
        all = all && one;
        any = any || one;
      }
      return (this->kind == Kind::all) ? all : any;
    }
    case Kind::acceptance:
    {
      double target;
      if (!adaptor.acceptance_target(target))
        Rcpp::stop("until = acceptance(…): the adaptor of " + adaptor.state_label() +
                   " has no target acceptance rate (am and rb_am do not; use stable(…) or iterations(…)).");
      this->acceptances.push_back(std::isnan(acceptance_probability) ? 0.0 : acceptance_probability);
      while (this->acceptances.size() > this->window)
        this->acceptances.pop_front();
      if (this->acceptances.size() < this->window)
        return false;
      double sum = 0.0;
      for (double a : this->acceptances)
        sum = sum + a;
      return std::abs(sum/double(this->window) - target) <= this->tolerance;
    }
    case Kind::stable:
    {
      std::vector<std::string> names = this->watch;
      if (names.empty() || ((names.size() == 1) && (names[0] == "default")))
        names = adaptor.default_watch();
      std::vector<double> now;
      for (const std::string &name : names)
      {
        std::vector<double> values;
        if (!adaptor.watched(name, values))
        {
          std::string available;
          for (const std::string &a : adaptor.watchable())
            available = available + (available.empty() ? "" : ", ") + a;
          Rcpp::stop("until = stable(…): the adaptor of " + adaptor.state_label() + " cannot watch " + name +
                     " (it can watch: " + available + ").");
        }
        now.insert(now.end(), values.begin(), values.end());
      }
      if (!this->history.empty() && (this->history.back().size() != now.size()))
        this->history.clear();
      this->history.push_back(now);
      while (this->history.size() > this->window + 1)
        this->history.pop_front();
      if (this->history.size() < this->window + 1)
        return false;
      for (const std::vector<double> &before : this->history)
        for (size_t i = 0; i < now.size(); ++i)
        {
          double change = std::abs(now[i] - before[i]);
          if (this->relative)
            change = change/std::max(std::abs(now[i]), 1e-300);
          if (!(change <= this->tolerance))
            return false; // (a NaN is not stable)
        }
      return true;
    }
  }
  return false;
}

}
