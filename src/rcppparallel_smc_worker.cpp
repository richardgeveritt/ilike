#include "rcppparallel_smc_worker.h"

namespace ilike
{
//Default constructor.
RcppParallelSMCWorker::RcppParallelSMCWorker()
{
  this->log_unnormalised_incremental_weights = NULL;
}

RcppParallelSMCWorker::RcppParallelSMCWorker(SMC* the_smc_in,
                                             size_t grain_size_in)
:SMCWorker(the_smc_in)
{
  this->log_unnormalised_incremental_weights = NULL;
  this->simulate_worker = SimulateWorker(this);
  this->conditional_simulate_worker = ConditionalSimulateWorker(this);
  this->move_worker = MoveWorker(this);
  this->weight_worker = WeightWorker(this);
  this->pf_initial_weight_worker = PFInitialWeightWorker(this);
  this->smcfixed_weight_worker = SMCFixedWeightWorker(this);
  this->smcadaptive_given_smcfixed_weight_worker = SMCAdaptiveGivenSMCFixedWeightWorker(this);
  this->smcadaptive_given_smcfixed_evaluate_target_worker = SMCAdaptiveGivenSMCFixedEvaluateTargetWorker(this);
  this->marginal_weight_worker = MarginalWeightWorker(this);
  this->generic_weight_worker = GenericWeightWorker(this);
  this->pf_weight_worker = PFWeightWorker(this);
  this->grain_size = grain_size_in;
}

//Copy constructor for the RcppParallelSMCWorker class.
RcppParallelSMCWorker::RcppParallelSMCWorker(const RcppParallelSMCWorker &another)
:SMCWorker(another)
{
  this->make_copy(another);
}

//Destructor for the RcppParallelSMCWorker class.
RcppParallelSMCWorker::~RcppParallelSMCWorker()
{
}

void RcppParallelSMCWorker::operator=(const RcppParallelSMCWorker &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  SMCWorker::operator=(another);
  this->make_copy(another);
}

SMCWorker* RcppParallelSMCWorker::duplicate() const
{
  return( new RcppParallelSMCWorker(*this));
}

void RcppParallelSMCWorker::make_copy(const RcppParallelSMCWorker &another)
{
  this->simulate_worker = another.simulate_worker;
  this->conditional_simulate_worker = another.conditional_simulate_worker;
  this->move_worker = another.move_worker;
  this->weight_worker = another.weight_worker;
  this->pf_initial_weight_worker = another.pf_initial_weight_worker;
  this->smcfixed_weight_worker = another.smcfixed_weight_worker;
  this->smcadaptive_given_smcfixed_weight_worker = another.smcadaptive_given_smcfixed_weight_worker;
  this->smcadaptive_given_smcfixed_evaluate_target_worker = another.smcadaptive_given_smcfixed_evaluate_target_worker;
  this->marginal_weight_worker = another.marginal_weight_worker;
  this->generic_weight_worker = another.generic_weight_worker;
  this->pf_weight_worker = another.pf_weight_worker;
  
  this->subsample_simulate_worker = another.subsample_simulate_worker;
  this->subsample_conditional_simulate_worker = another.subsample_conditional_simulate_worker;
  this->subsample_move_worker = another.subsample_move_worker;
  this->subsample_weight_worker = another.subsample_weight_worker;
  this->subsample_pf_initial_weight_worker = another.subsample_pf_initial_weight_worker;
  this->subsample_smcfixed_weight_worker = another.subsample_smcfixed_weight_worker;
  this->subsample_smcadaptive_given_smcfixed_weight_worker = another.subsample_smcadaptive_given_smcfixed_weight_worker;
  this->subsample_smcadaptive_given_smcfixed_evaluate_target_worker = another.subsample_smcadaptive_given_smcfixed_evaluate_target_worker;
  this->subsample_marginal_weight_worker = another.subsample_marginal_weight_worker;
  this->subsample_generic_weight_worker = another.subsample_generic_weight_worker;
  this->subsample_pf_weight_worker = another.subsample_pf_weight_worker;
  
  this->grain_size = another.grain_size;
  this->log_unnormalised_incremental_weights = another.log_unnormalised_incremental_weights;
}

arma::colvec RcppParallelSMCWorker::get_unnormalised_log_incremental_weights() const
{
  if (this->log_unnormalised_incremental_weights!=NULL)
    return arma::conv_to<arma::colvec>::from(*this->log_unnormalised_incremental_weights);
  else
    Rcpp::stop("RcppParallelSMCWorker::get_unnormalised_log_incremental_weights - weights not yet set.");
}

void RcppParallelSMCWorker::specific_simulate(Particles* next_particles)
{
  this->simulate_worker.particles_pointer = &next_particles->particles;
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    this->simulate_worker.particles_pointer->push_back(NULL);
  }
  parallelFor(0, this->get_number_of_particles(), this->simulate_worker, this->grain_size);
}

void RcppParallelSMCWorker::specific_simulate(Particles* next_particles,
                                              const Parameters &conditioned_on_parameters)
{
  this->conditional_simulate_worker.particles_pointer = &next_particles->particles;
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    this->conditional_simulate_worker.particles_pointer->push_back(NULL);
  }
  this->conditional_simulate_worker.conditioned_on_parameters_pointer = &conditioned_on_parameters;
  parallelFor(0, this->get_number_of_particles(), this->conditional_simulate_worker, this->grain_size);
}

void RcppParallelSMCWorker::subsample_specific_simulate(Particles* next_particles)
{
  this->subsample_simulate_worker.particles_pointer = &next_particles->particles;
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    this->subsample_simulate_worker.particles_pointer->push_back(NULL);
  }
  parallelFor(0, this->get_number_of_particles(), this->subsample_simulate_worker, this->grain_size);
}

void RcppParallelSMCWorker::subsample_specific_simulate(Particles* next_particles,
                                                        const Parameters &conditioned_on_parameters)
{
  this->subsample_conditional_simulate_worker.particles_pointer = &next_particles->particles;
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    this->subsample_conditional_simulate_worker.particles_pointer->push_back(NULL);
  }
  this->conditional_simulate_worker.conditioned_on_parameters_pointer = &conditioned_on_parameters;
  parallelFor(0, this->get_number_of_particles(), this->subsample_conditional_simulate_worker, this->grain_size);
}

void RcppParallelSMCWorker::specific_move(Particles* next_particles,
                                          const Particles* current_particles)
{
  this->move_worker.particles_pointer = &next_particles->particles;
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    this->move_worker.particles_pointer->push_back(NULL);
  }
  
  this->move_worker.current_particles_pointer = current_particles;
  parallelFor(0, this->get_number_of_particles(), this->move_worker, this->grain_size);
}

void RcppParallelSMCWorker::subsample_specific_move(Particles* next_particles,
                                                    const Particles* current_particles)
{
  this->subsample_move_worker.particles_pointer = &next_particles->particles;
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    this->subsample_move_worker.particles_pointer->push_back(NULL);
  }
  this->subsample_move_worker.current_particles_pointer = current_particles;
  parallelFor(0, this->get_number_of_particles(), this->subsample_move_worker, this->grain_size);
}

void RcppParallelSMCWorker::weight(const Index* index,
                                   Particles &current_particles)
{
  this->weight_worker.index_pointer = index;
  this->weight_worker.current_particles_pointer = &current_particles;
  this->weight_worker.my_log_unnormalised_incremental_weights = std::vector<double>(this->get_number_of_particles());
  
  parallelFor(0, this->get_number_of_particles(), this->weight_worker, this->grain_size);
  this->log_unnormalised_incremental_weights = &this->weight_worker.my_log_unnormalised_incremental_weights;
}

void RcppParallelSMCWorker::pf_initial_weight(Particles &current_particles)
{
  this->pf_initial_weight_worker.current_particles_pointer = &current_particles;
  this->pf_initial_weight_worker.my_log_unnormalised_incremental_weights = std::vector<double>(this->get_number_of_particles());
  
  parallelFor(0, this->get_number_of_particles(), this->pf_initial_weight_worker, this->grain_size);
  this->log_unnormalised_incremental_weights = &this->pf_initial_weight_worker.my_log_unnormalised_incremental_weights;
}

void RcppParallelSMCWorker::smcfixed_weight(const Index* index,
                                            Particles &current_particles)
{
  this->smcfixed_weight_worker.index_pointer = index;
  this->smcfixed_weight_worker.current_particles_pointer = &current_particles.particles;
  this->smcfixed_weight_worker.my_log_unnormalised_incremental_weights = std::vector<double>(this->get_number_of_particles());
  
  parallelFor(0, this->get_number_of_particles(), this->smcfixed_weight_worker, this->grain_size);
  this->log_unnormalised_incremental_weights = &this->smcfixed_weight_worker.my_log_unnormalised_incremental_weights;
}

void RcppParallelSMCWorker::smcadaptive_given_smcfixed_weight(const Index* index,
                                                              Particles &current_particles)
{
  this->smcadaptive_given_smcfixed_weight_worker.index_pointer = index;
  this->smcadaptive_given_smcfixed_weight_worker.current_particles_pointer = &current_particles.particles;
  this->smcadaptive_given_smcfixed_weight_worker.my_log_unnormalised_incremental_weights = std::vector<double>(this->get_number_of_particles());
  
  parallelFor(0, this->get_number_of_particles(), this->smcadaptive_given_smcfixed_weight_worker, this->grain_size);
  this->log_unnormalised_incremental_weights = &this->smcadaptive_given_smcfixed_weight_worker.my_log_unnormalised_incremental_weights;
}

void RcppParallelSMCWorker::smcadaptive_given_smcfixed_evaluate_target(const Index* index,
                                                                       Particles &current_particles)
{
  this->smcadaptive_given_smcfixed_evaluate_target_worker.index_pointer = index;
  this->smcadaptive_given_smcfixed_evaluate_target_worker.current_particles_pointer = &current_particles;
  this->smcadaptive_given_smcfixed_weight_worker.my_log_unnormalised_incremental_weights = std::vector<double>(this->get_number_of_particles());
  
  parallelFor(0, this->get_number_of_particles(), this->smcadaptive_given_smcfixed_evaluate_target_worker, this->grain_size);
  this->log_unnormalised_incremental_weights = &this->smcadaptive_given_smcfixed_weight_worker.my_log_unnormalised_incremental_weights;
}

void RcppParallelSMCWorker::marginal_weight(const Index* index,
                                            Particles &current_particles,
                                            Particles &previous_particles,
                                            ProposalKernel* proposal_kernel)
{
  Rcpp::stop("Don't use yet - not thread safe.");
  this->marginal_weight_worker.index_pointer = index;
  this->marginal_weight_worker.current_particles_pointer = &current_particles;
  this->marginal_weight_worker.previous_particles_pointer = &previous_particles;
  this->marginal_weight_worker.proposal_kernel_pointer = proposal_kernel;
  this->marginal_weight_worker.my_log_unnormalised_incremental_weights = std::vector<double>(this->get_number_of_particles());
  
  parallelFor(0, this->get_number_of_particles(), this->marginal_weight_worker, this->grain_size);
  this->log_unnormalised_incremental_weights = &this->marginal_weight_worker.my_log_unnormalised_incremental_weights;
}

void RcppParallelSMCWorker::generic_weight(const Index* index,
                                           Particles &current_particles,
                                           Particles &previous_particles,
                                           ProposalKernel* proposal_kernel,
                                           ProposalKernel* L_kernel)
{
  this->generic_weight_worker.index_pointer = index;
  this->generic_weight_worker.current_particles_pointer = &current_particles;
  this->generic_weight_worker.previous_particles_pointer = &previous_particles;
  this->generic_weight_worker.proposal_kernel_pointer = proposal_kernel;
  this->generic_weight_worker.L_kernel_pointer = L_kernel;
  this->generic_weight_worker.my_log_unnormalised_incremental_weights = std::vector<double>(this->get_number_of_particles());;
  
  parallelFor(0, this->get_number_of_particles(), this->generic_weight_worker, this->grain_size);
  this->log_unnormalised_incremental_weights = &this->generic_weight_worker.my_log_unnormalised_incremental_weights;
}

void RcppParallelSMCWorker::pf_weight(const Index* index,
                                      Particles &current_particles,
                                      Particles &previous_particles,
                                      ProposalKernel* proposal_kernel)
{
  this->pf_weight_worker.index_pointer = index;
  this->pf_weight_worker.current_particles_pointer = &current_particles;
  this->pf_weight_worker.previous_particles_pointer = &previous_particles;
  this->pf_weight_worker.proposal_kernel_pointer = proposal_kernel;
  this->pf_weight_worker.my_log_unnormalised_incremental_weights = std::vector<double>(this->get_number_of_particles());;
  
  parallelFor(0, this->get_number_of_particles(), this->pf_weight_worker, this->grain_size);
  this->log_unnormalised_incremental_weights = &this->pf_weight_worker.my_log_unnormalised_incremental_weights;
}

void RcppParallelSMCWorker::subsample_weight(const Index* index,
                                             Particles &current_particles)
{
  this->subsample_weight_worker.index_pointer = index;
  this->subsample_weight_worker.current_particles_pointer = &current_particles;
  this->subsample_weight_worker.my_log_unnormalised_incremental_weights = std::vector<double>(this->get_number_of_particles());;
  
  parallelFor(0, this->get_number_of_particles(), this->subsample_weight_worker, this->grain_size);
  this->log_unnormalised_incremental_weights = &this->subsample_weight_worker.my_log_unnormalised_incremental_weights;
}

void RcppParallelSMCWorker::subsample_pf_initial_weight(Particles &current_particles)
{
  this->subsample_pf_initial_weight_worker.current_particles_pointer = &current_particles;
  this->subsample_pf_weight_worker.my_log_unnormalised_incremental_weights = std::vector<double>(this->get_number_of_particles());;
  
  parallelFor(0, this->get_number_of_particles(), this->subsample_pf_initial_weight_worker, this->grain_size);
  this->log_unnormalised_incremental_weights = &this->subsample_pf_weight_worker.my_log_unnormalised_incremental_weights;
}

void RcppParallelSMCWorker::subsample_smcfixed_weight(const Index* index,
                                                      Particles &current_particles)
{
  this->subsample_smcfixed_weight_worker.index_pointer = index;
  this->subsample_smcfixed_weight_worker.current_particles_pointer = &current_particles;
  this->subsample_smcfixed_weight_worker.my_log_unnormalised_incremental_weights = std::vector<double>(this->get_number_of_particles());;
  
  parallelFor(0, this->get_number_of_particles(), this->subsample_smcfixed_weight_worker, this->grain_size);
  this->log_unnormalised_incremental_weights = &this->subsample_smcfixed_weight_worker.my_log_unnormalised_incremental_weights;
}

void RcppParallelSMCWorker::subsample_smcadaptive_given_smcfixed_weight(const Index* index,
                                                                        Particles &current_particles)
{
  this->subsample_smcadaptive_given_smcfixed_weight_worker.index_pointer = index;
  this->subsample_smcadaptive_given_smcfixed_weight_worker.current_particles_pointer = &current_particles;
  this->subsample_smcadaptive_given_smcfixed_weight_worker.my_log_unnormalised_incremental_weights = std::vector<double>(this->get_number_of_particles());;
  
  parallelFor(0, this->get_number_of_particles(), this->subsample_smcadaptive_given_smcfixed_weight_worker, this->grain_size);
  this->log_unnormalised_incremental_weights = &this->subsample_smcadaptive_given_smcfixed_weight_worker.my_log_unnormalised_incremental_weights;
}

void RcppParallelSMCWorker::subsample_smcadaptive_given_smcfixed_evaluate_target(const Index* index,
                                                                                 Particles &current_particles)
{
  this->subsample_smcadaptive_given_smcfixed_evaluate_target_worker.index_pointer = index;
  this->subsample_smcadaptive_given_smcfixed_evaluate_target_worker.current_particles_pointer = &current_particles;
  this->subsample_smcadaptive_given_smcfixed_weight_worker.my_log_unnormalised_incremental_weights = std::vector<double>(this->get_number_of_particles());;
  
  parallelFor(0, this->get_number_of_particles(), this->subsample_smcadaptive_given_smcfixed_evaluate_target_worker, this->grain_size);
  this->log_unnormalised_incremental_weights = &this->subsample_smcadaptive_given_smcfixed_weight_worker.my_log_unnormalised_incremental_weights;
}

void RcppParallelSMCWorker::subsample_marginal_weight(const Index* index,
                                                      Particles &current_particles,
                                                      Particles &previous_particles,
                                                      ProposalKernel* proposal_kernel)
{
  this->subsample_marginal_weight_worker.index_pointer = index;
  this->subsample_marginal_weight_worker.current_particles_pointer = &current_particles;
  this->subsample_marginal_weight_worker.previous_particles_pointer = &previous_particles;
  this->subsample_marginal_weight_worker.proposal_kernel_pointer = proposal_kernel;
  this->subsample_marginal_weight_worker.my_log_unnormalised_incremental_weights = std::vector<double>(this->get_number_of_particles());;
  
  parallelFor(0, this->get_number_of_particles(), this->subsample_marginal_weight_worker, this->grain_size);
  this->log_unnormalised_incremental_weights = &this->subsample_marginal_weight_worker.my_log_unnormalised_incremental_weights;
}

void RcppParallelSMCWorker::subsample_generic_weight(const Index* index,
                                                     Particles &current_particles,
                                                     Particles &previous_particles,
                                                     ProposalKernel* proposal_kernel,
                                                     ProposalKernel* L_kernel)
{
  this->subsample_generic_weight_worker.index_pointer = index;
  this->subsample_generic_weight_worker.current_particles_pointer = &current_particles;
  this->subsample_generic_weight_worker.previous_particles_pointer = &previous_particles;
  this->subsample_generic_weight_worker.proposal_kernel_pointer = proposal_kernel;
  this->subsample_generic_weight_worker.L_kernel_pointer = L_kernel;
  this->subsample_generic_weight_worker.my_log_unnormalised_incremental_weights = std::vector<double>(this->get_number_of_particles());;
  
  parallelFor(0, this->get_number_of_particles(), this->subsample_generic_weight_worker, this->grain_size);
  this->log_unnormalised_incremental_weights = &this->subsample_generic_weight_worker.my_log_unnormalised_incremental_weights;
}

void RcppParallelSMCWorker::subsample_pf_weight(const Index* index,
                                                Particles &current_particles,
                                                Particles &previous_particles,
                                                ProposalKernel* proposal_kernel)
{
  this->subsample_pf_weight_worker.index_pointer = index;
  this->subsample_pf_weight_worker.current_particles_pointer = &current_particles;
  this->subsample_pf_weight_worker.previous_particles_pointer = &previous_particles;
  this->subsample_pf_weight_worker.proposal_kernel_pointer = proposal_kernel;
  this->subsample_pf_weight_worker.my_log_unnormalised_incremental_weights = std::vector<double>(this->get_number_of_particles());;
  
  parallelFor(0, this->get_number_of_particles(), this->subsample_pf_weight_worker, this->grain_size);
  this->log_unnormalised_incremental_weights = &this->subsample_pf_weight_worker.my_log_unnormalised_incremental_weights;
}
}
