#include "kalman_updater.h"
#include "kalman_filter_output.h"
#include "kalman_filter.h"

KalmanUpdater::KalmanUpdater()
{
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
}
