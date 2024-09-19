#include "kalman_updater.h"
#include "kalman_filter_output.h"
#include "kalman_filter.h"

namespace ilike
{
KalmanUpdater::KalmanUpdater()
{
  
}

KalmanUpdater::KalmanUpdater(const std::string &state_variable_in,
                             const std::string &measurement_variable_in)
{
  this->state_variable = state_variable_in;
  this->measurement_variable = measurement_variable_in;
}

KalmanUpdater::~KalmanUpdater()
{
  
}

KalmanUpdater::KalmanUpdater(const KalmanUpdater &another)
{
  this->make_copy(another);
}

void KalmanUpdater::operator=(const KalmanUpdater &another)
{
  if(this == &another)
    return;
  
  this->make_copy(another);
}

void KalmanUpdater::make_copy(const KalmanUpdater &another)
{
  this->set_using_parameters = another.set_using_parameters;
  this->conditioned_on_parameters = another.conditioned_on_parameters;
  
  this->state_variable = another.state_variable;
  this->measurement_variable = another.measurement_variable;
}

std::string KalmanUpdater::get_state_variable() const
{
  return this->state_variable;
}

std::string KalmanUpdater::get_measurement_variable() const
{
  return this->measurement_variable;
}
}
