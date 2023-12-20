#include "kalman_predictor.h"
#include "kalman_filter_output.h"

KalmanPredictor::KalmanPredictor()
{
}

KalmanPredictor::~KalmanPredictor()
{

}

KalmanPredictor::KalmanPredictor(const KalmanPredictor &another)
{
  this->make_copy(another);
}

void KalmanPredictor::operator=(const KalmanPredictor &another)
{
  if(this == &another)
    return;

  this->make_copy(another);
}

void KalmanPredictor::make_copy(const KalmanPredictor &another)
{
  this->set_using_parameters = another.set_using_parameters;
  //this->set_using_time = another.set_using_time;
  this->conditioned_on_parameters = another.conditioned_on_parameters;
}
