#include "rcppparallel_workers.h"
#include "particle_simulator.h"
#include "smc.h"
#include "single_point_move_output.h"
#include "rcppparallel_smc_worker.h"
#include "vector_single_index.h"
#include "utils.h"

SimulateWorker::SimulateWorker()
{
}

SimulateWorker::SimulateWorker(RcppParallelSMCWorker* smc_worker_in)
{
  this->smc_worker = smc_worker_in;
  this->particles_pointer = NULL;
}

SimulateWorker::~SimulateWorker()
{

}

SimulateWorker::SimulateWorker(const SimulateWorker &another)
{
  this->make_copy(another);
}

void SimulateWorker::operator=(const SimulateWorker &another)
{
  if(this == &another){ //if a==a
    return;
  }

  this->make_copy(another);
}

void SimulateWorker::operator()(std::size_t begin, std::size_t end)
{
  RandomNumberGenerator local_rng(*this->smc_worker->get_rng());
  local_rng.seed(this->smc_worker->get_seed(),end);
  for (std::size_t i = begin; i < end; ++i)
  {
    (*this->particles_pointer)[i] = new SinglePointMoveOutput(this->smc_worker->the_smc->particle_simulator->simulate(local_rng,
                                                                              this->smc_worker->the_smc->factors,
                                                                              &this->smc_worker->the_smc->proposals_to_transform_for,
                                                                              &this->smc_worker->the_smc->proposals_to_find_gradient_for,
                                                                              this->smc_worker->the_smc->sequencer.schedule_parameters));

    if (this->smc_worker->the_smc->proposal_is_evaluated==true)
    {
      (*this->particles_pointer)[i]->back().previous_target_evaluated = this->smc_worker->the_smc->particle_simulator->evaluate((*this->particles_pointer)[i]->back());
    }
  }
}

void SimulateWorker::make_copy(const SimulateWorker &another)
{
  this->smc_worker = another.smc_worker;
  this->particles_pointer = another.particles_pointer;
}


ConditionalSimulateWorker::ConditionalSimulateWorker()
{
}

ConditionalSimulateWorker::ConditionalSimulateWorker(RcppParallelSMCWorker* smc_worker_in)
{
  this->smc_worker = smc_worker_in;
  this->particles_pointer = NULL;
  this->conditioned_on_parameters_pointer = NULL;
}

ConditionalSimulateWorker::~ConditionalSimulateWorker()
{

}

ConditionalSimulateWorker::ConditionalSimulateWorker(const ConditionalSimulateWorker &another)
{
  this->make_copy(another);
}

void ConditionalSimulateWorker::operator=(const ConditionalSimulateWorker &another)
{
  if(this == &another){ //if a==a
    return;
  }

  this->make_copy(another);
}

void ConditionalSimulateWorker::operator()(std::size_t begin, std::size_t end)
{
  RandomNumberGenerator local_rng(*this->smc_worker->get_rng());
  local_rng.seed(this->smc_worker->get_seed(),end);
  for (std::size_t i = begin; i < end; ++i)
  {
    (*this->particles_pointer)[i] = new SinglePointMoveOutput(this->smc_worker->the_smc->particle_simulator->simulate(local_rng,
                                                                                                                      this->smc_worker->the_smc->factors,
                                                                                                                      &this->smc_worker->the_smc->proposals_to_transform_for,
                                                                                                                      &this->smc_worker->the_smc->proposals_to_find_gradient_for,
                                                                                                                      *this->conditioned_on_parameters_pointer,
                                                                                                                      this->smc_worker->the_smc->sequencer.schedule_parameters));

    if (this->smc_worker->the_smc->proposal_is_evaluated==true)
    {
      (*this->particles_pointer)[i]->back().previous_target_evaluated = this->smc_worker->the_smc->particle_simulator->evaluate((*this->particles_pointer)[i]->back());
    }
  }
}

void ConditionalSimulateWorker::make_copy(const ConditionalSimulateWorker &another)
{
  this->smc_worker = another.smc_worker;
  this->particles_pointer = another.particles_pointer;
  this->conditioned_on_parameters_pointer = another.conditioned_on_parameters_pointer;
}


SubsampleSimulateWorker::SubsampleSimulateWorker()
{
}

SubsampleSimulateWorker::SubsampleSimulateWorker(RcppParallelSMCWorker* smc_worker_in)
{
  this->smc_worker = smc_worker_in;
  this->particles_pointer = NULL;
}

SubsampleSimulateWorker::~SubsampleSimulateWorker()
{

}

SubsampleSimulateWorker::SubsampleSimulateWorker(const SubsampleSimulateWorker &another)
{
  this->make_copy(another);
}

void SubsampleSimulateWorker::operator=(const SubsampleSimulateWorker &another)
{
  if(this == &another){ //if a==a
    return;
  }

  this->make_copy(another);
}

void SubsampleSimulateWorker::operator()(std::size_t begin, std::size_t end)
{
  RandomNumberGenerator local_rng(*this->smc_worker->get_rng());
  local_rng.seed(this->smc_worker->get_seed(),end);
  for (std::size_t i = begin; i < end; ++i)
  {
    (*this->particles_pointer)[i] = new SinglePointMoveOutput(this->smc_worker->the_smc->particle_simulator->subsample_simulate(local_rng,
                                                                                                                      this->smc_worker->the_smc->factors,
                                                                                                                      &this->smc_worker->the_smc->proposals_to_transform_for,
                                                                                                                      &this->smc_worker->the_smc->proposals_to_find_gradient_for,
                                                                                                                      this->smc_worker->the_smc->sequencer.schedule_parameters));

    if (this->smc_worker->the_smc->proposal_is_evaluated==true)
    {
      (*this->particles_pointer)[i]->back().previous_target_evaluated = this->smc_worker->the_smc->particle_simulator->subsample_evaluate((*this->particles_pointer)[i]->back());
    }
  }
}

void SubsampleSimulateWorker::make_copy(const SubsampleSimulateWorker &another)
{
  this->smc_worker = another.smc_worker;
  this->particles_pointer = another.particles_pointer;
}


SubsampleConditionalSimulateWorker::SubsampleConditionalSimulateWorker()
{
}

SubsampleConditionalSimulateWorker::SubsampleConditionalSimulateWorker(RcppParallelSMCWorker* smc_worker_in)
{
  this->smc_worker = smc_worker_in;
  this->particles_pointer = NULL;
  this->conditioned_on_parameters_pointer = NULL;
}

SubsampleConditionalSimulateWorker::~SubsampleConditionalSimulateWorker()
{

}

SubsampleConditionalSimulateWorker::SubsampleConditionalSimulateWorker(const SubsampleConditionalSimulateWorker &another)
{
  this->make_copy(another);
}

void SubsampleConditionalSimulateWorker::operator=(const SubsampleConditionalSimulateWorker &another)
{
  if(this == &another){ //if a==a
    return;
  }

  this->make_copy(another);
}

void SubsampleConditionalSimulateWorker::operator()(std::size_t begin, std::size_t end)
{
  RandomNumberGenerator local_rng(*this->smc_worker->get_rng());
  local_rng.seed(this->smc_worker->get_seed(),end);
  for (std::size_t i = begin; i < end; ++i)
  {
    (*this->particles_pointer)[i] = new SinglePointMoveOutput(this->smc_worker->the_smc->particle_simulator->subsample_simulate(local_rng,
                                                                                                                      this->smc_worker->the_smc->factors,
                                                                                                                      &this->smc_worker->the_smc->proposals_to_transform_for,
                                                                                                                      &this->smc_worker->the_smc->proposals_to_find_gradient_for,
                                                                                                                      *this->conditioned_on_parameters_pointer,
                                                                                                                      this->smc_worker->the_smc->sequencer.schedule_parameters));

    if (this->smc_worker->the_smc->proposal_is_evaluated==true)
    {
      (*this->particles_pointer)[i]->back().previous_target_evaluated = this->smc_worker->the_smc->particle_simulator->subsample_evaluate((*this->particles_pointer)[i]->back());
    }
  }
}

void SubsampleConditionalSimulateWorker::make_copy(const SubsampleConditionalSimulateWorker &another)
{
  this->smc_worker = another.smc_worker;
  this->particles_pointer = another.particles_pointer;
  this->conditioned_on_parameters_pointer = another.conditioned_on_parameters_pointer;
}


MoveWorker::MoveWorker()
{
}

MoveWorker::MoveWorker(RcppParallelSMCWorker* smc_worker_in)
{
  this->smc_worker = smc_worker_in;
  this->particles_pointer = NULL;
  this->current_particles_pointer = NULL;
}

MoveWorker::~MoveWorker()
{

}

MoveWorker::MoveWorker(const MoveWorker &another)
{
  this->make_copy(another);
}

void MoveWorker::operator=(const MoveWorker &another)
{
  if(this == &another){ //if a==a
    return;
  }

  this->make_copy(another);
}

void MoveWorker::operator()(std::size_t begin, std::size_t end)
{
  RandomNumberGenerator local_rng(*this->smc_worker->get_rng());
  local_rng.seed(this->smc_worker->get_seed(),end);

  for (std::size_t i = begin; i < end; ++i)
  {
    (*this->particles_pointer)[i] = this->smc_worker->the_smc->move(local_rng,
                                                                    (*this->current_particles_pointer)[current_particles_pointer->ancestor_variables[i]]->back());
  }
}

void MoveWorker::make_copy(const MoveWorker &another)
{
  this->smc_worker = another.smc_worker;
  this->particles_pointer = another.particles_pointer;
  this->current_particles_pointer = another.current_particles_pointer;
}


SubsampleMoveWorker::SubsampleMoveWorker()
{
}

SubsampleMoveWorker::SubsampleMoveWorker(RcppParallelSMCWorker* smc_worker_in)
{
  this->smc_worker = smc_worker_in;
  this->particles_pointer = NULL;
  this->current_particles_pointer = NULL;
}

SubsampleMoveWorker::~SubsampleMoveWorker()
{

}

SubsampleMoveWorker::SubsampleMoveWorker(const SubsampleMoveWorker &another)
{
  this->make_copy(another);
}

void SubsampleMoveWorker::operator=(const SubsampleMoveWorker &another)
{
  if(this == &another){ //if a==a
    return;
  }

  this->make_copy(another);
}

void SubsampleMoveWorker::operator()(std::size_t begin, std::size_t end)
{
  RandomNumberGenerator local_rng(*this->smc_worker->get_rng());
  local_rng.seed(this->smc_worker->get_seed(),end);
  for (std::size_t i = begin; i < end; ++i)
  {
    (*this->particles_pointer)[i] = this->smc_worker->the_smc->subsample_move(local_rng,
                                                                              (*this->current_particles_pointer)[current_particles_pointer->ancestor_variables[i]]->back());
  }
}

void SubsampleMoveWorker::make_copy(const SubsampleMoveWorker &another)
{
  this->smc_worker = another.smc_worker;
  this->particles_pointer = another.particles_pointer;
  this->current_particles_pointer = another.current_particles_pointer;
}


WeightWorker::WeightWorker()
{
}

WeightWorker::WeightWorker(RcppParallelSMCWorker* smc_worker_in)
{
  this->smc_worker = smc_worker_in;
  this->current_particles_pointer = NULL;
  this->index_pointer = NULL;
}

WeightWorker::~WeightWorker()
{

}

WeightWorker::WeightWorker(const WeightWorker &another)
{
  this->make_copy(another);
}

void WeightWorker::operator=(const WeightWorker &another)
{
  if(this == &another){ //if a==a
    return;
  }

  this->make_copy(another);
}

void WeightWorker::operator()(std::size_t begin, std::size_t end)
{
  for (std::size_t i = begin; i < end; ++i)
  {
    this->my_log_unnormalised_incremental_weights[i] = (*this->current_particles_pointer)[i]->back().evaluate_likelihoods(this->index_pointer) - (*this->current_particles_pointer)[i]->back().previous_target_evaluated;
  }
}

void WeightWorker::make_copy(const WeightWorker &another)
{
  this->smc_worker = another.smc_worker;
  this->my_log_unnormalised_incremental_weights = another.my_log_unnormalised_incremental_weights;
  this->index_pointer = another.index_pointer;
  this->current_particles_pointer = another.current_particles_pointer;
}


PFInitialWeightWorker::PFInitialWeightWorker()
{
}

PFInitialWeightWorker::PFInitialWeightWorker(RcppParallelSMCWorker* smc_worker_in)
{
  this->smc_worker = smc_worker_in;
  this->current_particles_pointer = NULL;
}

PFInitialWeightWorker::~PFInitialWeightWorker()
{

}

PFInitialWeightWorker::PFInitialWeightWorker(const PFInitialWeightWorker &another)
{
  this->make_copy(another);
}

void PFInitialWeightWorker::operator=(const PFInitialWeightWorker &another)
{
  if(this == &another){ //if a==a
    return;
  }

  this->make_copy(another);
}

void PFInitialWeightWorker::operator()(std::size_t begin, std::size_t end)
{
  VectorSingleIndex index(0);

  for (std::size_t i = begin; i < end; ++i)
  {
    this->my_log_unnormalised_incremental_weights[i] = (*this->current_particles_pointer)[i]->back().evaluate_likelihoods(&index) - (*this->current_particles_pointer)[i]->back().previous_target_evaluated;
  }
}

void PFInitialWeightWorker::make_copy(const PFInitialWeightWorker &another)
{
  this->smc_worker = another.smc_worker;
  this->my_log_unnormalised_incremental_weights = another.my_log_unnormalised_incremental_weights;
  this->current_particles_pointer = another.current_particles_pointer;
}


SMCFixedWeightWorker::SMCFixedWeightWorker()
{
}

SMCFixedWeightWorker::SMCFixedWeightWorker(RcppParallelSMCWorker* smc_worker_in)
{
  this->smc_worker = smc_worker_in;
  this->current_particles_pointer = NULL;
  this->index_pointer = NULL;
}

SMCFixedWeightWorker::~SMCFixedWeightWorker()
{

}

SMCFixedWeightWorker::SMCFixedWeightWorker(const SMCFixedWeightWorker &another)
{
  this->make_copy(another);
}

void SMCFixedWeightWorker::operator=(const SMCFixedWeightWorker &another)
{
  if(this == &another){ //if a==a
    return;
  }

  this->make_copy(another);
}

void SMCFixedWeightWorker::operator()(std::size_t begin, std::size_t end)
{
  for (std::size_t i = begin; i < end; ++i)
  {
    (*this->current_particles_pointer)[i]->back().evaluate_smcfixed_part_of_likelihoods(this->index_pointer);
  }
}

void SMCFixedWeightWorker::make_copy(const SMCFixedWeightWorker &another)
{
  this->smc_worker = another.smc_worker;
  this->my_log_unnormalised_incremental_weights = another.my_log_unnormalised_incremental_weights;
  this->index_pointer = another.index_pointer;
  this->current_particles_pointer = another.current_particles_pointer;
}


SMCAdaptiveGivenSMCFixedWeightWorker::SMCAdaptiveGivenSMCFixedWeightWorker()
{
}

SMCAdaptiveGivenSMCFixedWeightWorker::SMCAdaptiveGivenSMCFixedWeightWorker(RcppParallelSMCWorker* smc_worker_in)
{
  this->smc_worker = smc_worker_in;
  this->current_particles_pointer = NULL;
  this->index_pointer = NULL;
}

SMCAdaptiveGivenSMCFixedWeightWorker::~SMCAdaptiveGivenSMCFixedWeightWorker()
{

}

SMCAdaptiveGivenSMCFixedWeightWorker::SMCAdaptiveGivenSMCFixedWeightWorker(const SMCAdaptiveGivenSMCFixedWeightWorker &another)
{
  this->make_copy(another);
}

void SMCAdaptiveGivenSMCFixedWeightWorker::operator=(const SMCAdaptiveGivenSMCFixedWeightWorker &another)
{
  if(this == &another){ //if a==a
    return;
  }

  this->make_copy(another);
}

void SMCAdaptiveGivenSMCFixedWeightWorker::operator()(std::size_t begin, std::size_t end)
{
  for (std::size_t i = begin; i < end; ++i)
  {
    //double prev = current_particles[i]->back().previous_target_evaluated;
    double a = (*this->current_particles_pointer)[i]->back().evaluate_smcadaptive_part_given_smcfixed_likelihoods(this->index_pointer);
    double b = (*this->current_particles_pointer)[i]->back().previous_target_evaluated;
    if (a==-arma::datum::inf)
    {
      this->my_log_unnormalised_incremental_weights[i] = -arma::datum::inf;
    }
    else
    {
      this->my_log_unnormalised_incremental_weights[i] = a - b;
    }
    
    //this->my_log_unnormalised_incremental_weights[i] = (*this->current_particles_pointer)[i]->back().evaluate_smcadaptive_part_given_smcfixed_likelihoods(this->index_pointer) - (*this->current_particles_pointer)[i]->back().previous_target_evaluated;
  }
}

void SMCAdaptiveGivenSMCFixedWeightWorker::make_copy(const SMCAdaptiveGivenSMCFixedWeightWorker &another)
{
  this->smc_worker = another.smc_worker;
  this->my_log_unnormalised_incremental_weights = another.my_log_unnormalised_incremental_weights;
  this->index_pointer = another.index_pointer;
  this->current_particles_pointer = another.current_particles_pointer;
}


SMCAdaptiveGivenSMCFixedEvaluateTargetWorker::SMCAdaptiveGivenSMCFixedEvaluateTargetWorker()
{
}

SMCAdaptiveGivenSMCFixedEvaluateTargetWorker::SMCAdaptiveGivenSMCFixedEvaluateTargetWorker(RcppParallelSMCWorker* smc_worker_in)
{
  this->smc_worker = smc_worker_in;
  this->current_particles_pointer = NULL;
  this->index_pointer = NULL;
}

SMCAdaptiveGivenSMCFixedEvaluateTargetWorker::~SMCAdaptiveGivenSMCFixedEvaluateTargetWorker()
{

}

SMCAdaptiveGivenSMCFixedEvaluateTargetWorker::SMCAdaptiveGivenSMCFixedEvaluateTargetWorker(const SMCAdaptiveGivenSMCFixedEvaluateTargetWorker &another)
{
  this->make_copy(another);
}

void SMCAdaptiveGivenSMCFixedEvaluateTargetWorker::operator=(const SMCAdaptiveGivenSMCFixedEvaluateTargetWorker &another)
{
  if(this == &another){ //if a==a
    return;
  }

  this->make_copy(another);
}

void SMCAdaptiveGivenSMCFixedEvaluateTargetWorker::operator()(std::size_t begin, std::size_t end)
{
  for (std::size_t i = begin; i < end; ++i)
  {
    this->my_log_unnormalised_incremental_weights[i] = (*this->current_particles_pointer)[i]->back().evaluate_smcadaptive_part_given_smcfixed_likelihoods(this->index_pointer);
  }
}

void SMCAdaptiveGivenSMCFixedEvaluateTargetWorker::make_copy(const SMCAdaptiveGivenSMCFixedEvaluateTargetWorker &another)
{
  this->smc_worker = another.smc_worker;
  this->my_log_unnormalised_incremental_weights = another.my_log_unnormalised_incremental_weights;
  this->index_pointer = another.index_pointer;
  this->current_particles_pointer = another.current_particles_pointer;
}


MarginalWeightWorker::MarginalWeightWorker()
{
}

MarginalWeightWorker::MarginalWeightWorker(RcppParallelSMCWorker* smc_worker_in)
{
  this->smc_worker = smc_worker_in;
  this->current_particles_pointer = NULL;
  this->previous_particles_pointer = NULL;
  this->proposal_kernel_pointer = NULL;
  this->index_pointer = NULL;
}

MarginalWeightWorker::~MarginalWeightWorker()
{

}

MarginalWeightWorker::MarginalWeightWorker(const MarginalWeightWorker &another)
{
  this->make_copy(another);
}

void MarginalWeightWorker::operator=(const MarginalWeightWorker &another)
{
  if(this == &another){ //if a==a
    return;
  }

  this->make_copy(another);
}

void MarginalWeightWorker::operator()(std::size_t begin, std::size_t end)
{

  for (std::size_t i = begin; i < end; ++i)
  {
    arma::colvec terms(this->smc_worker->get_number_of_particles());

    // If proposal is Gaussian, might be able to use the fast Gauss transform in low-dim
    for (size_t j = 0; j < this->smc_worker->get_number_of_particles(); ++j)
    {
      terms[j] = this->current_particles_pointer->previous_normalised_log_weights[j] + this->proposal_kernel_pointer->evaluate_kernel((*this->current_particles_pointer)[i]->back(),
                                                                                                         (*this->previous_particles_pointer)[this->previous_particles_pointer->ancestor_variables[j]]->back());
    }

    this->my_log_unnormalised_incremental_weights[i] = (*this->current_particles_pointer)[i]->back().evaluate_likelihoods(this->index_pointer) - log_sum_exp(terms);
  }
}

void MarginalWeightWorker::make_copy(const MarginalWeightWorker &another)
{
  this->smc_worker = another.smc_worker;
  this->my_log_unnormalised_incremental_weights = another.my_log_unnormalised_incremental_weights;
  this->index_pointer = another.index_pointer;
  this->current_particles_pointer = another.current_particles_pointer;
  this->previous_particles_pointer = another.previous_particles_pointer;
  this->proposal_kernel_pointer = another.proposal_kernel_pointer;
}


GenericWeightWorker::GenericWeightWorker()
{
}

GenericWeightWorker::GenericWeightWorker(RcppParallelSMCWorker* smc_worker_in)
{
  this->smc_worker = smc_worker_in;
  this->current_particles_pointer = NULL;
  this->previous_particles_pointer = NULL;
  this->proposal_kernel_pointer = NULL;
  this->L_kernel_pointer = NULL;
  this->index_pointer = NULL;
}

GenericWeightWorker::~GenericWeightWorker()
{

}

GenericWeightWorker::GenericWeightWorker(const GenericWeightWorker &another)
{
  this->make_copy(another);
}

void GenericWeightWorker::operator=(const GenericWeightWorker &another)
{
  if(this == &another){ //if a==a
    return;
  }

  this->make_copy(another);
}

void GenericWeightWorker::operator()(std::size_t begin, std::size_t end)
{
  for (std::size_t i = begin; i < end; ++i)
  {
    this->my_log_unnormalised_incremental_weights[i] = (*this->current_particles_pointer)[i]->back().evaluate_likelihoods(this->index_pointer)
    + this->L_kernel_pointer->evaluate_kernel((*this->previous_particles_pointer)[this->previous_particles_pointer->ancestor_variables[i]]->back(),
                                (*this->current_particles_pointer)[i]->back())
    - (*this->previous_particles_pointer)[this->previous_particles_pointer->ancestor_variables[i]]->back().target_evaluated
    - this->proposal_kernel_pointer->evaluate_kernel((*this->current_particles_pointer)[i]->back(),
                                       (*this->previous_particles_pointer)[this->previous_particles_pointer->ancestor_variables[i]]->back());
  }
}

void GenericWeightWorker::make_copy(const GenericWeightWorker &another)
{
  this->smc_worker = another.smc_worker;
  this->my_log_unnormalised_incremental_weights = another.my_log_unnormalised_incremental_weights;
  this->index_pointer = another.index_pointer;
  this->current_particles_pointer = another.current_particles_pointer;
  this->previous_particles_pointer = another.previous_particles_pointer;
  this->proposal_kernel_pointer = another.proposal_kernel_pointer;
  this->L_kernel_pointer = another.L_kernel_pointer;
}


PFWeightWorker::PFWeightWorker()
{
}

PFWeightWorker::PFWeightWorker(RcppParallelSMCWorker* smc_worker_in)
{
  this->smc_worker = smc_worker_in;
  this->current_particles_pointer = NULL;
  this->previous_particles_pointer = NULL;
  this->proposal_kernel_pointer = NULL;
  this->index_pointer = NULL;
}

PFWeightWorker::~PFWeightWorker()
{

}

PFWeightWorker::PFWeightWorker(const PFWeightWorker &another)
{
  this->make_copy(another);
}

void PFWeightWorker::operator=(const PFWeightWorker &another)
{
  if(this == &another){ //if a==a
    return;
  }

  this->make_copy(another);
}

void PFWeightWorker::operator()(std::size_t begin, std::size_t end)
{
  for (std::size_t i = begin; i < end; ++i)
  {
    this->my_log_unnormalised_incremental_weights[i] = (*this->current_particles_pointer)[i]->back().evaluate_likelihoods(this->index_pointer);
    if (this->proposal_kernel_pointer!=NULL)
      this->my_log_unnormalised_incremental_weights[i] = this->my_log_unnormalised_incremental_weights[i] - this->proposal_kernel_pointer->evaluate_kernel((*this->current_particles_pointer)[i]->back(),
                                                                                                                                                                           (*this->previous_particles_pointer)[this->previous_particles_pointer->ancestor_variables[i]]->back());
  }
}

void PFWeightWorker::make_copy(const PFWeightWorker &another)
{
  this->smc_worker = another.smc_worker;
  this->my_log_unnormalised_incremental_weights = another.my_log_unnormalised_incremental_weights;
  this->index_pointer = another.index_pointer;
  this->current_particles_pointer = another.current_particles_pointer;
  this->previous_particles_pointer = another.previous_particles_pointer;
  this->proposal_kernel_pointer = another.proposal_kernel_pointer;
}


SubsampleWeightWorker::SubsampleWeightWorker()
{
}

SubsampleWeightWorker::SubsampleWeightWorker(RcppParallelSMCWorker* smc_worker_in)
{
  this->smc_worker = smc_worker_in;
  this->current_particles_pointer = NULL;
  this->index_pointer = NULL;
}

SubsampleWeightWorker::~SubsampleWeightWorker()
{

}

SubsampleWeightWorker::SubsampleWeightWorker(const SubsampleWeightWorker &another)
{
  this->make_copy(another);
}

void SubsampleWeightWorker::operator=(const SubsampleWeightWorker &another)
{
  if(this == &another){ //if a==a
    return;
  }

  this->make_copy(another);
}

void SubsampleWeightWorker::operator()(std::size_t begin, std::size_t end)
{
  for (std::size_t i = begin; i < end; ++i)
  {
    this->my_log_unnormalised_incremental_weights[i] = (*this->current_particles_pointer)[i]->back().subsample_evaluate_likelihoods(this->index_pointer) - (*this->current_particles_pointer)[i]->back().subsample_previous_target_evaluated;
  }
}

void SubsampleWeightWorker::make_copy(const SubsampleWeightWorker &another)
{
  this->smc_worker = another.smc_worker;
  this->my_log_unnormalised_incremental_weights = another.my_log_unnormalised_incremental_weights;
  this->index_pointer = another.index_pointer;
  this->current_particles_pointer = another.current_particles_pointer;
}


SubsamplePFInitialWeightWorker::SubsamplePFInitialWeightWorker()
{
}

SubsamplePFInitialWeightWorker::SubsamplePFInitialWeightWorker(RcppParallelSMCWorker* smc_worker_in)
{
  this->smc_worker = smc_worker_in;
  this->current_particles_pointer = NULL;
}

SubsamplePFInitialWeightWorker::~SubsamplePFInitialWeightWorker()
{

}

SubsamplePFInitialWeightWorker::SubsamplePFInitialWeightWorker(const SubsamplePFInitialWeightWorker &another)
{
  this->make_copy(another);
}

void SubsamplePFInitialWeightWorker::operator=(const SubsamplePFInitialWeightWorker &another)
{
  if(this == &another){ //if a==a
    return;
  }

  this->make_copy(another);
}

void SubsamplePFInitialWeightWorker::operator()(std::size_t begin, std::size_t end)
{
  VectorSingleIndex index(0);

  for (std::size_t i = begin; i < end; ++i)
  {
    this->my_log_unnormalised_incremental_weights[i] = (*this->current_particles_pointer)[i]->back().subsample_evaluate_likelihoods(&index) - (*this->current_particles_pointer)[i]->back().subsample_previous_target_evaluated;
  }
}

void SubsamplePFInitialWeightWorker::make_copy(const SubsamplePFInitialWeightWorker &another)
{
  this->smc_worker = another.smc_worker;
  this->my_log_unnormalised_incremental_weights = another.my_log_unnormalised_incremental_weights;
  this->current_particles_pointer = another.current_particles_pointer;
}


SubsampleSMCFixedWeightWorker::SubsampleSMCFixedWeightWorker()
{
}

SubsampleSMCFixedWeightWorker::SubsampleSMCFixedWeightWorker(RcppParallelSMCWorker* smc_worker_in)
{
  this->smc_worker = smc_worker_in;
  this->current_particles_pointer = NULL;
  this->index_pointer = NULL;
}

SubsampleSMCFixedWeightWorker::~SubsampleSMCFixedWeightWorker()
{

}

SubsampleSMCFixedWeightWorker::SubsampleSMCFixedWeightWorker(const SubsampleSMCFixedWeightWorker &another)
{
  this->make_copy(another);
}

void SubsampleSMCFixedWeightWorker::operator=(const SubsampleSMCFixedWeightWorker &another)
{
  if(this == &another){ //if a==a
    return;
  }

  this->make_copy(another);
}

void SubsampleSMCFixedWeightWorker::operator()(std::size_t begin, std::size_t end)
{
  for (std::size_t i = begin; i < end; ++i)
  {
    (*this->current_particles_pointer)[i]->back().subsample_evaluate_smcfixed_part_of_likelihoods(this->index_pointer);
  }
}

void SubsampleSMCFixedWeightWorker::make_copy(const SubsampleSMCFixedWeightWorker &another)
{
  this->smc_worker = another.smc_worker;
  this->my_log_unnormalised_incremental_weights = another.my_log_unnormalised_incremental_weights;
  this->index_pointer = another.index_pointer;
  this->current_particles_pointer = another.current_particles_pointer;
}


SubsampleSMCAdaptiveGivenSMCFixedWeightWorker::SubsampleSMCAdaptiveGivenSMCFixedWeightWorker()
{
}

SubsampleSMCAdaptiveGivenSMCFixedWeightWorker::SubsampleSMCAdaptiveGivenSMCFixedWeightWorker(RcppParallelSMCWorker* smc_worker_in)
{
  this->smc_worker = smc_worker_in;
  this->current_particles_pointer = NULL;
  this->index_pointer = NULL;
}

SubsampleSMCAdaptiveGivenSMCFixedWeightWorker::~SubsampleSMCAdaptiveGivenSMCFixedWeightWorker()
{

}

SubsampleSMCAdaptiveGivenSMCFixedWeightWorker::SubsampleSMCAdaptiveGivenSMCFixedWeightWorker(const SubsampleSMCAdaptiveGivenSMCFixedWeightWorker &another)
{
  this->make_copy(another);
}

void SubsampleSMCAdaptiveGivenSMCFixedWeightWorker::operator=(const SubsampleSMCAdaptiveGivenSMCFixedWeightWorker &another)
{
  if(this == &another){ //if a==a
    return;
  }

  this->make_copy(another);
}

void SubsampleSMCAdaptiveGivenSMCFixedWeightWorker::operator()(std::size_t begin, std::size_t end)
{
  for (std::size_t i = begin; i < end; ++i)
  {
    this->my_log_unnormalised_incremental_weights[i] = (*this->current_particles_pointer)[i]->back().subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(this->index_pointer) - (*this->current_particles_pointer)[i]->back().subsample_previous_target_evaluated;
  }
}

void SubsampleSMCAdaptiveGivenSMCFixedWeightWorker::make_copy(const SubsampleSMCAdaptiveGivenSMCFixedWeightWorker &another)
{
  this->smc_worker = another.smc_worker;
  this->my_log_unnormalised_incremental_weights = another.my_log_unnormalised_incremental_weights;
  this->index_pointer = another.index_pointer;
  this->current_particles_pointer = another.current_particles_pointer;
}


SubsampleSMCAdaptiveGivenSMCFixedEvaluateTargetWorker::SubsampleSMCAdaptiveGivenSMCFixedEvaluateTargetWorker()
{
}

SubsampleSMCAdaptiveGivenSMCFixedEvaluateTargetWorker::SubsampleSMCAdaptiveGivenSMCFixedEvaluateTargetWorker(RcppParallelSMCWorker* smc_worker_in)
{
  this->smc_worker = smc_worker_in;
  this->current_particles_pointer = NULL;
  this->index_pointer = NULL;
}

SubsampleSMCAdaptiveGivenSMCFixedEvaluateTargetWorker::~SubsampleSMCAdaptiveGivenSMCFixedEvaluateTargetWorker()
{

}

SubsampleSMCAdaptiveGivenSMCFixedEvaluateTargetWorker::SubsampleSMCAdaptiveGivenSMCFixedEvaluateTargetWorker(const SubsampleSMCAdaptiveGivenSMCFixedEvaluateTargetWorker &another)
{
  this->make_copy(another);
}

void SubsampleSMCAdaptiveGivenSMCFixedEvaluateTargetWorker::operator=(const SubsampleSMCAdaptiveGivenSMCFixedEvaluateTargetWorker &another)
{
  if(this == &another){ //if a==a
    return;
  }

  this->make_copy(another);
}

void SubsampleSMCAdaptiveGivenSMCFixedEvaluateTargetWorker::operator()(std::size_t begin, std::size_t end)
{
  for (std::size_t i = begin; i < end; ++i)
  {
    this->my_log_unnormalised_incremental_weights[i] = (*this->current_particles_pointer)[i]->back().subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(this->index_pointer);
  }
}

void SubsampleSMCAdaptiveGivenSMCFixedEvaluateTargetWorker::make_copy(const SubsampleSMCAdaptiveGivenSMCFixedEvaluateTargetWorker &another)
{
  this->smc_worker = another.smc_worker;
  this->my_log_unnormalised_incremental_weights = another.my_log_unnormalised_incremental_weights;
  this->index_pointer = another.index_pointer;
  this->current_particles_pointer = another.current_particles_pointer;
}


SubsampleMarginalWeightWorker::SubsampleMarginalWeightWorker()
{
}

SubsampleMarginalWeightWorker::SubsampleMarginalWeightWorker(RcppParallelSMCWorker* smc_worker_in)
{
  this->smc_worker = smc_worker_in;
  this->current_particles_pointer = NULL;
  this->previous_particles_pointer = NULL;
  this->proposal_kernel_pointer = NULL;
  this->index_pointer = NULL;
}

SubsampleMarginalWeightWorker::~SubsampleMarginalWeightWorker()
{

}

SubsampleMarginalWeightWorker::SubsampleMarginalWeightWorker(const SubsampleMarginalWeightWorker &another)
{
  this->make_copy(another);
}

void SubsampleMarginalWeightWorker::operator=(const SubsampleMarginalWeightWorker &another)
{
  if(this == &another){ //if a==a
    return;
  }

  this->make_copy(another);
}

void SubsampleMarginalWeightWorker::operator()(std::size_t begin, std::size_t end)
{

  for (std::size_t i = begin; i < end; ++i)
  {
    arma::colvec terms(this->smc_worker->get_number_of_particles());

    // If proposal is Gaussian, might be able to use the fast Gauss transform in low-dim
    for (size_t j = 0; j < this->smc_worker->get_number_of_particles(); ++j)
    {
      terms[j] = this->current_particles_pointer->previous_normalised_log_weights[j] + this->proposal_kernel_pointer->subsample_evaluate_kernel((*this->current_particles_pointer)[i]->back(),
                                                                                                                                      (*this->previous_particles_pointer)[this->previous_particles_pointer->ancestor_variables[j]]->back());
    }

    this->my_log_unnormalised_incremental_weights[i] = (*this->current_particles_pointer)[i]->back().subsample_evaluate_likelihoods(this->index_pointer) - log_sum_exp(terms);
  }
}

void SubsampleMarginalWeightWorker::make_copy(const SubsampleMarginalWeightWorker &another)
{
  this->smc_worker = another.smc_worker;
  this->my_log_unnormalised_incremental_weights = another.my_log_unnormalised_incremental_weights;
  this->index_pointer = another.index_pointer;
  this->current_particles_pointer = another.current_particles_pointer;
  this->previous_particles_pointer = another.previous_particles_pointer;
  this->proposal_kernel_pointer = another.proposal_kernel_pointer;
}


SubsampleGenericWeightWorker::SubsampleGenericWeightWorker()
{
}

SubsampleGenericWeightWorker::SubsampleGenericWeightWorker(RcppParallelSMCWorker* smc_worker_in)
{
  this->smc_worker = smc_worker_in;
  this->current_particles_pointer = NULL;
  this->previous_particles_pointer = NULL;
  this->proposal_kernel_pointer = NULL;
  this->L_kernel_pointer = NULL;
  this->index_pointer = NULL;
}

SubsampleGenericWeightWorker::~SubsampleGenericWeightWorker()
{

}

SubsampleGenericWeightWorker::SubsampleGenericWeightWorker(const SubsampleGenericWeightWorker &another)
{
  this->make_copy(another);
}

void SubsampleGenericWeightWorker::operator=(const SubsampleGenericWeightWorker &another)
{
  if(this == &another){ //if a==a
    return;
  }

  this->make_copy(another);
}

void SubsampleGenericWeightWorker::operator()(std::size_t begin, std::size_t end)
{
  for (std::size_t i = begin; i < end; ++i)
  {
    this->my_log_unnormalised_incremental_weights[i] = (*this->current_particles_pointer)[i]->back().subsample_evaluate_likelihoods(this->index_pointer)
    + this->L_kernel_pointer->subsample_evaluate_kernel((*this->previous_particles_pointer)[this->previous_particles_pointer->ancestor_variables[i]]->back(),
                                              (*this->current_particles_pointer)[i]->back())
    - (*this->previous_particles_pointer)[this->previous_particles_pointer->ancestor_variables[i]]->back().subsample_target_evaluated
    - this->proposal_kernel_pointer->subsample_evaluate_kernel((*this->current_particles_pointer)[i]->back(),
                                                     (*this->previous_particles_pointer)[this->previous_particles_pointer->ancestor_variables[i]]->back());
  }
}

void SubsampleGenericWeightWorker::make_copy(const SubsampleGenericWeightWorker &another)
{
  this->smc_worker = another.smc_worker;
  this->my_log_unnormalised_incremental_weights = another.my_log_unnormalised_incremental_weights;
  this->index_pointer = another.index_pointer;
  this->current_particles_pointer = another.current_particles_pointer;
  this->previous_particles_pointer = another.previous_particles_pointer;
  this->proposal_kernel_pointer = another.proposal_kernel_pointer;
  this->L_kernel_pointer = another.L_kernel_pointer;
}


SubsamplePFWeightWorker::SubsamplePFWeightWorker()
{
}

SubsamplePFWeightWorker::SubsamplePFWeightWorker(RcppParallelSMCWorker* smc_worker_in)
{
  this->smc_worker = smc_worker_in;
  this->current_particles_pointer = NULL;
  this->previous_particles_pointer = NULL;
  this->proposal_kernel_pointer = NULL;
  this->index_pointer = NULL;
}

SubsamplePFWeightWorker::~SubsamplePFWeightWorker()
{

}

SubsamplePFWeightWorker::SubsamplePFWeightWorker(const SubsamplePFWeightWorker &another)
{
  this->make_copy(another);
}

void SubsamplePFWeightWorker::operator=(const SubsamplePFWeightWorker &another)
{
  if(this == &another){ //if a==a
    return;
  }

  this->make_copy(another);
}

void SubsamplePFWeightWorker::operator()(std::size_t begin, std::size_t end)
{
  for (std::size_t i = begin; i < end; ++i)
  {
    this->my_log_unnormalised_incremental_weights[i] = (*this->current_particles_pointer)[i]->back().subsample_evaluate_likelihoods(this->index_pointer);
    if (this->proposal_kernel_pointer!=NULL)
      this->my_log_unnormalised_incremental_weights[i] = this->my_log_unnormalised_incremental_weights[i] - this->proposal_kernel_pointer->subsample_evaluate_kernel((*this->current_particles_pointer)[i]->back(),
                                                                                                                                                                           (*this->previous_particles_pointer)[this->previous_particles_pointer->ancestor_variables[i]]->back());
  }
}

void SubsamplePFWeightWorker::make_copy(const SubsamplePFWeightWorker &another)
{
  this->smc_worker = another.smc_worker;
  this->my_log_unnormalised_incremental_weights = another.my_log_unnormalised_incremental_weights;
  this->index_pointer = another.index_pointer;
  this->current_particles_pointer = another.current_particles_pointer;
  this->previous_particles_pointer = another.previous_particles_pointer;
  this->proposal_kernel_pointer = another.proposal_kernel_pointer;
}
