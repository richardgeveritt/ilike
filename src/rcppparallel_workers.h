#ifndef RCPPPARALLELWORKERS_H
#define RCPPPARALLELWORKERS_H

#include <RcppParallel.h>

#include <vector>

#include "particles.h"

class MoveOutput;

class RcppParallelSMCWorker;

class SimulateWorker : public RcppParallel::Worker {

public:

  SimulateWorker();

  SimulateWorker(RcppParallelSMCWorker* smc_worker_in);

  ~SimulateWorker();

  SimulateWorker(const SimulateWorker &another);

  void operator=(const SimulateWorker &another);

  void operator()(std::size_t begin, std::size_t end);

  std::vector< MoveOutput* >* particles_pointer;

private:

  friend class RcppParallelSMCWorker;
  RcppParallelSMCWorker* smc_worker;

  void make_copy(const SimulateWorker &another);

};

class ConditionalSimulateWorker : public RcppParallel::Worker {
  
public:
  
  ConditionalSimulateWorker();
  
  ConditionalSimulateWorker(RcppParallelSMCWorker* smc_worker_in);
  
  ~ConditionalSimulateWorker();
  
  ConditionalSimulateWorker(const ConditionalSimulateWorker &another);
  
  void operator=(const ConditionalSimulateWorker &another);
  
  void operator()(std::size_t begin, std::size_t end);
  
  std::vector< MoveOutput* >* particles_pointer;
  
  const Parameters* conditioned_on_parameters_pointer;
  
private:
  
  friend class RcppParallelSMCWorker;
  RcppParallelSMCWorker* smc_worker;
  
  void make_copy(const ConditionalSimulateWorker &another);
  
};

class SubsampleSimulateWorker : public RcppParallel::Worker {
  
public:
  
  SubsampleSimulateWorker();
  
  SubsampleSimulateWorker(RcppParallelSMCWorker* smc_worker_in);
  
  ~SubsampleSimulateWorker();
  
  SubsampleSimulateWorker(const SubsampleSimulateWorker &another);
  
  void operator=(const SubsampleSimulateWorker &another);
  
  void operator()(std::size_t begin, std::size_t end);
  
  std::vector< MoveOutput* >* particles_pointer;
  
private:
  
  friend class RcppParallelSMCWorker;
  RcppParallelSMCWorker* smc_worker;
  
  void make_copy(const SubsampleSimulateWorker &another);
  
};

class SubsampleConditionalSimulateWorker : public RcppParallel::Worker {
  
public:
  
  SubsampleConditionalSimulateWorker();
  
  SubsampleConditionalSimulateWorker(RcppParallelSMCWorker* smc_worker_in);
  
  ~SubsampleConditionalSimulateWorker();
  
  SubsampleConditionalSimulateWorker(const SubsampleConditionalSimulateWorker &another);
  
  void operator=(const SubsampleConditionalSimulateWorker &another);
  
  void operator()(std::size_t begin, std::size_t end);
  
  std::vector< MoveOutput* >* particles_pointer;
  
  const Parameters* conditioned_on_parameters_pointer;
  
private:
  
  friend class RcppParallelSMCWorker;
  RcppParallelSMCWorker* smc_worker;
  
  void make_copy(const SubsampleConditionalSimulateWorker &another);
  
};


class MoveWorker : public RcppParallel::Worker {
  
public:
  
  MoveWorker();
  
  MoveWorker(RcppParallelSMCWorker* smc_worker_in);
  
  ~MoveWorker();
  
  MoveWorker(const MoveWorker &another);
  
  void operator=(const MoveWorker &another);
  
  void operator()(std::size_t begin, std::size_t end);
  
  std::vector< MoveOutput* >* particles_pointer;
  const Particles* current_particles_pointer;
  
private:
  
  friend class RcppParallelSMCWorker;
  RcppParallelSMCWorker* smc_worker;
  
  void make_copy(const MoveWorker &another);
  
};


class SubsampleMoveWorker : public RcppParallel::Worker {
  
public:
  
  SubsampleMoveWorker();
  
  SubsampleMoveWorker(RcppParallelSMCWorker* smc_worker_in);
  
  ~SubsampleMoveWorker();
  
  SubsampleMoveWorker(const SubsampleMoveWorker &another);
  
  void operator=(const SubsampleMoveWorker &another);
  
  void operator()(std::size_t begin, std::size_t end);
  
  std::vector< MoveOutput* >* particles_pointer;
  const Particles* current_particles_pointer;
  
private:
  
  friend class RcppParallelSMCWorker;
  RcppParallelSMCWorker* smc_worker;
  
  void make_copy(const SubsampleMoveWorker &another);
  
};


class WeightWorker : public RcppParallel::Worker {
  
public:
  
  WeightWorker();
  
  WeightWorker(RcppParallelSMCWorker* smc_worker_in);
  
  ~WeightWorker();
  
  WeightWorker(const WeightWorker &another);
  
  void operator=(const WeightWorker &another);
  
  void operator()(std::size_t begin, std::size_t end);
  
  std::vector<double> my_log_unnormalised_incremental_weights;
  const Index* index_pointer;
  Particles* current_particles_pointer;
  
private:
  
  friend class RcppParallelSMCWorker;
  RcppParallelSMCWorker* smc_worker;
  
  void make_copy(const WeightWorker &another);
  
};


class PFInitialWeightWorker : public RcppParallel::Worker {
  
public:
  
  PFInitialWeightWorker();
  
  PFInitialWeightWorker(RcppParallelSMCWorker* smc_worker_in);
  
  ~PFInitialWeightWorker();
  
  PFInitialWeightWorker(const PFInitialWeightWorker &another);
  
  void operator=(const PFInitialWeightWorker &another);
  
  void operator()(std::size_t begin, std::size_t end);
  
  std::vector<double> my_log_unnormalised_incremental_weights;
  Particles* current_particles_pointer;
  
private:
  
  friend class RcppParallelSMCWorker;
  RcppParallelSMCWorker* smc_worker;
  
  void make_copy(const PFInitialWeightWorker &another);
  
};


class SMCFixedWeightWorker : public RcppParallel::Worker {
  
public:
  
  SMCFixedWeightWorker();
  
  SMCFixedWeightWorker(RcppParallelSMCWorker* smc_worker_in);
  
  ~SMCFixedWeightWorker();
  
  SMCFixedWeightWorker(const SMCFixedWeightWorker &another);
  
  void operator=(const SMCFixedWeightWorker &another);
  
  void operator()(std::size_t begin, std::size_t end);
  
  std::vector<double> my_log_unnormalised_incremental_weights;
  const Index* index_pointer;
  std::vector< MoveOutput* >* current_particles_pointer;
  //Particles* current_particles_pointer;
  
private:
  
  friend class RcppParallelSMCWorker;
  RcppParallelSMCWorker* smc_worker;
  
  void make_copy(const SMCFixedWeightWorker &another);
  
};


class SMCAdaptiveGivenSMCFixedWeightWorker : public RcppParallel::Worker {
  
public:
  
  SMCAdaptiveGivenSMCFixedWeightWorker();
  
  SMCAdaptiveGivenSMCFixedWeightWorker(RcppParallelSMCWorker* smc_worker_in);
  
  ~SMCAdaptiveGivenSMCFixedWeightWorker();
  
  SMCAdaptiveGivenSMCFixedWeightWorker(const SMCAdaptiveGivenSMCFixedWeightWorker &another);
  
  void operator=(const SMCAdaptiveGivenSMCFixedWeightWorker &another);
  
  void operator()(std::size_t begin, std::size_t end);
  
  std::vector<double> my_log_unnormalised_incremental_weights;
  const Index* index_pointer;
  std::vector< MoveOutput* >* current_particles_pointer;
  
private:
  
  friend class RcppParallelSMCWorker;
  RcppParallelSMCWorker* smc_worker;
  
  void make_copy(const SMCAdaptiveGivenSMCFixedWeightWorker &another);
  
};


class SMCAdaptiveGivenSMCFixedEvaluateTargetWorker : public RcppParallel::Worker {
  
public:
  
  SMCAdaptiveGivenSMCFixedEvaluateTargetWorker();
  
  SMCAdaptiveGivenSMCFixedEvaluateTargetWorker(RcppParallelSMCWorker* smc_worker_in);
  
  ~SMCAdaptiveGivenSMCFixedEvaluateTargetWorker();
  
  SMCAdaptiveGivenSMCFixedEvaluateTargetWorker(const SMCAdaptiveGivenSMCFixedEvaluateTargetWorker &another);
  
  void operator=(const SMCAdaptiveGivenSMCFixedEvaluateTargetWorker &another);
  
  void operator()(std::size_t begin, std::size_t end);
  
  std::vector<double> my_log_unnormalised_incremental_weights;
  const Index* index_pointer;
  Particles* current_particles_pointer;
  
private:
  
  friend class RcppParallelSMCWorker;
  RcppParallelSMCWorker* smc_worker;
  
  void make_copy(const SMCAdaptiveGivenSMCFixedEvaluateTargetWorker &another);
  
};


class MarginalWeightWorker : public RcppParallel::Worker {
  
public:
  
  MarginalWeightWorker();
  
  MarginalWeightWorker(RcppParallelSMCWorker* smc_worker_in);
  
  ~MarginalWeightWorker();
  
  MarginalWeightWorker(const MarginalWeightWorker &another);
  
  void operator=(const MarginalWeightWorker &another);
  
  void operator()(std::size_t begin, std::size_t end);
  
  std::vector<double> my_log_unnormalised_incremental_weights;
  const Index* index_pointer;
  Particles* current_particles_pointer;
  Particles* previous_particles_pointer;
  ProposalKernel* proposal_kernel_pointer;
  
private:
  
  friend class RcppParallelSMCWorker;
  RcppParallelSMCWorker* smc_worker;
  
  void make_copy(const MarginalWeightWorker &another);
  
};


class GenericWeightWorker : public RcppParallel::Worker {
  
public:
  
  GenericWeightWorker();
  
  GenericWeightWorker(RcppParallelSMCWorker* smc_worker_in);
  
  ~GenericWeightWorker();
  
  GenericWeightWorker(const GenericWeightWorker &another);
  
  void operator=(const GenericWeightWorker &another);
  
  void operator()(std::size_t begin, std::size_t end);
  
  std::vector<double> my_log_unnormalised_incremental_weights;
  const Index* index_pointer;
  Particles* current_particles_pointer;
  Particles* previous_particles_pointer;
  ProposalKernel* proposal_kernel_pointer;
  ProposalKernel* L_kernel_pointer;
  
private:
  
  friend class RcppParallelSMCWorker;
  RcppParallelSMCWorker* smc_worker;
  
  void make_copy(const GenericWeightWorker &another);
  
};

class PFWeightWorker : public RcppParallel::Worker {
  
public:
  
  PFWeightWorker();
  
  PFWeightWorker(RcppParallelSMCWorker* smc_worker_in);
  
  ~PFWeightWorker();
  
  PFWeightWorker(const PFWeightWorker &another);
  
  void operator=(const PFWeightWorker &another);
  
  void operator()(std::size_t begin, std::size_t end);
  
  std::vector<double> my_log_unnormalised_incremental_weights;
  const Index* index_pointer;
  Particles* current_particles_pointer;
  Particles* previous_particles_pointer;
  ProposalKernel* proposal_kernel_pointer;
  
private:
  
  friend class RcppParallelSMCWorker;
  RcppParallelSMCWorker* smc_worker;
  
  void make_copy(const PFWeightWorker &another);
  
};

class SubsampleWeightWorker : public RcppParallel::Worker {
  
public:
  
  SubsampleWeightWorker();
  
  SubsampleWeightWorker(RcppParallelSMCWorker* smc_worker_in);
  
  ~SubsampleWeightWorker();
  
  SubsampleWeightWorker(const SubsampleWeightWorker &another);
  
  void operator=(const SubsampleWeightWorker &another);
  
  void operator()(std::size_t begin, std::size_t end);
  
  std::vector<double> my_log_unnormalised_incremental_weights;
  const Index* index_pointer;
  Particles* current_particles_pointer;
  
private:
  
  friend class RcppParallelSMCWorker;
  RcppParallelSMCWorker* smc_worker;
  
  void make_copy(const SubsampleWeightWorker &another);
  
};


class SubsamplePFInitialWeightWorker : public RcppParallel::Worker {
  
public:
  
  SubsamplePFInitialWeightWorker();
  
  SubsamplePFInitialWeightWorker(RcppParallelSMCWorker* smc_worker_in);
  
  ~SubsamplePFInitialWeightWorker();
  
  SubsamplePFInitialWeightWorker(const SubsamplePFInitialWeightWorker &another);
  
  void operator=(const SubsamplePFInitialWeightWorker &another);
  
  void operator()(std::size_t begin, std::size_t end);
  
  std::vector<double> my_log_unnormalised_incremental_weights;
  Particles* current_particles_pointer;
  
private:
  
  friend class RcppParallelSMCWorker;
  RcppParallelSMCWorker* smc_worker;
  
  void make_copy(const SubsamplePFInitialWeightWorker &another);
  
};


class SubsampleSMCFixedWeightWorker : public RcppParallel::Worker {
  
public:
  
  SubsampleSMCFixedWeightWorker();
  
  SubsampleSMCFixedWeightWorker(RcppParallelSMCWorker* smc_worker_in);
  
  ~SubsampleSMCFixedWeightWorker();
  
  SubsampleSMCFixedWeightWorker(const SubsampleSMCFixedWeightWorker &another);
  
  void operator=(const SubsampleSMCFixedWeightWorker &another);
  
  void operator()(std::size_t begin, std::size_t end);
  
  std::vector<double> my_log_unnormalised_incremental_weights;
  const Index* index_pointer;
  Particles* current_particles_pointer;
  
private:
  
  friend class RcppParallelSMCWorker;
  RcppParallelSMCWorker* smc_worker;
  
  void make_copy(const SubsampleSMCFixedWeightWorker &another);
  
};


class SubsampleSMCAdaptiveGivenSMCFixedWeightWorker : public RcppParallel::Worker {
  
public:
  
  SubsampleSMCAdaptiveGivenSMCFixedWeightWorker();
  
  SubsampleSMCAdaptiveGivenSMCFixedWeightWorker(RcppParallelSMCWorker* smc_worker_in);
  
  ~SubsampleSMCAdaptiveGivenSMCFixedWeightWorker();
  
  SubsampleSMCAdaptiveGivenSMCFixedWeightWorker(const SubsampleSMCAdaptiveGivenSMCFixedWeightWorker &another);
  
  void operator=(const SubsampleSMCAdaptiveGivenSMCFixedWeightWorker &another);
  
  void operator()(std::size_t begin, std::size_t end);
  
  std::vector<double> my_log_unnormalised_incremental_weights;
  const Index* index_pointer;
  Particles* current_particles_pointer;
  
private:
  
  friend class RcppParallelSMCWorker;
  RcppParallelSMCWorker* smc_worker;
  
  void make_copy(const SubsampleSMCAdaptiveGivenSMCFixedWeightWorker &another);
  
};


class SubsampleSMCAdaptiveGivenSMCFixedEvaluateTargetWorker : public RcppParallel::Worker {
  
public:
  
  SubsampleSMCAdaptiveGivenSMCFixedEvaluateTargetWorker();
  
  SubsampleSMCAdaptiveGivenSMCFixedEvaluateTargetWorker(RcppParallelSMCWorker* smc_worker_in);
  
  ~SubsampleSMCAdaptiveGivenSMCFixedEvaluateTargetWorker();
  
  SubsampleSMCAdaptiveGivenSMCFixedEvaluateTargetWorker(const SubsampleSMCAdaptiveGivenSMCFixedEvaluateTargetWorker &another);
  
  void operator=(const SubsampleSMCAdaptiveGivenSMCFixedEvaluateTargetWorker &another);
  
  void operator()(std::size_t begin, std::size_t end);
  
  std::vector<double> my_log_unnormalised_incremental_weights;
  const Index* index_pointer;
  Particles* current_particles_pointer;
  
private:
  
  friend class RcppParallelSMCWorker;
  RcppParallelSMCWorker* smc_worker;
  
  void make_copy(const SubsampleSMCAdaptiveGivenSMCFixedEvaluateTargetWorker &another);
  
};


class SubsampleMarginalWeightWorker : public RcppParallel::Worker {
  
public:
  
  SubsampleMarginalWeightWorker();
  
  SubsampleMarginalWeightWorker(RcppParallelSMCWorker* smc_worker_in);
  
  ~SubsampleMarginalWeightWorker();
  
  SubsampleMarginalWeightWorker(const SubsampleMarginalWeightWorker &another);
  
  void operator=(const SubsampleMarginalWeightWorker &another);
  
  void operator()(std::size_t begin, std::size_t end);
  
  std::vector<double> my_log_unnormalised_incremental_weights;
  const Index* index_pointer;
  Particles* current_particles_pointer;
  Particles* previous_particles_pointer;
  ProposalKernel* proposal_kernel_pointer;
  
private:
  
  friend class RcppParallelSMCWorker;
  RcppParallelSMCWorker* smc_worker;
  
  void make_copy(const SubsampleMarginalWeightWorker &another);
  
};


class SubsampleGenericWeightWorker : public RcppParallel::Worker {
  
public:
  
  SubsampleGenericWeightWorker();
  
  SubsampleGenericWeightWorker(RcppParallelSMCWorker* smc_worker_in);
  
  ~SubsampleGenericWeightWorker();
  
  SubsampleGenericWeightWorker(const SubsampleGenericWeightWorker &another);
  
  void operator=(const SubsampleGenericWeightWorker &another);
  
  void operator()(std::size_t begin, std::size_t end);
  
  std::vector<double> my_log_unnormalised_incremental_weights;
  const Index* index_pointer;
  Particles* current_particles_pointer;
  Particles* previous_particles_pointer;
  ProposalKernel* proposal_kernel_pointer;
  ProposalKernel* L_kernel_pointer;
  
private:
  
  friend class RcppParallelSMCWorker;
  RcppParallelSMCWorker* smc_worker;
  
  void make_copy(const SubsampleGenericWeightWorker &another);
  
};

class SubsamplePFWeightWorker : public RcppParallel::Worker {
  
public:
  
  SubsamplePFWeightWorker();
  
  SubsamplePFWeightWorker(RcppParallelSMCWorker* smc_worker_in);
  
  ~SubsamplePFWeightWorker();
  
  SubsamplePFWeightWorker(const SubsamplePFWeightWorker &another);
  
  void operator=(const SubsamplePFWeightWorker &another);
  
  void operator()(std::size_t begin, std::size_t end);
  
  std::vector<double> my_log_unnormalised_incremental_weights;
  const Index* index_pointer;
  Particles* current_particles_pointer;
  Particles* previous_particles_pointer;
  ProposalKernel* proposal_kernel_pointer;
  
private:
  
  friend class RcppParallelSMCWorker;
  RcppParallelSMCWorker* smc_worker;
  
  void make_copy(const SubsamplePFWeightWorker &another);
  
};

#endif
