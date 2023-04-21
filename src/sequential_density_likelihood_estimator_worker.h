#ifndef SEQUENTIALDENSITYLIKELIHOODESTIMATORWORKER_H
#define SEQUENTIALDENSITYLIKELIHOODESTIMATORWORKER_H

#include "density_likelihood_estimator_worker.h"

#include <vector>
#include "parameters.h"

class DensityLikelihoodEstimatorOutput;

class SequentialDensityLikelihoodEstimatorWorker : public DensityLikelihoodEstimatorWorker
{
public:

  SequentialDensityLikelihoodEstimatorWorker(void);
  virtual ~SequentialDensityLikelihoodEstimatorWorker(void);

  SequentialDensityLikelihoodEstimatorWorker(DensityLikelihoodEstimator* the_dle_in);

  SequentialDensityLikelihoodEstimatorWorker(const SequentialDensityLikelihoodEstimatorWorker &another);
  void operator=(const SequentialDensityLikelihoodEstimatorWorker &another);
  DensityLikelihoodEstimatorWorker* duplicate() const;
  
  std::vector<Parameters> get_points() const;

protected:

  void specific_simulate(DensityLikelihoodEstimatorOutput* output);
  void specific_simulate(DensityLikelihoodEstimatorOutput* output,
                         const Parameters &conditioned_on_parameters);
  
  void subsample_specific_simulate(DensityLikelihoodEstimatorOutput* output);
  void subsample_specific_simulate(DensityLikelihoodEstimatorOutput* output,
                                   const Parameters &conditioned_on_parameters);

  void make_copy(const SequentialDensityLikelihoodEstimatorWorker &another);

  std::vector<Parameters> points;

};

#endif
