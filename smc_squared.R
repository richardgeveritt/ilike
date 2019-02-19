library(LaplacesDemon)
library(rapportools)
library(gsubfn)
library(svMisc)
library(ENmisc)
require(parallel)
library(lqmm)

logsumexp <- function(x){
  xmax = which.max(x)
  if (sum(xmax)==0)
  {
    -Inf
  }
  else
  {
    result = log1p(sum(exp(x[-xmax]-x[xmax])))+x[xmax]
    if (is.nan(result))
    {
      -Inf
    }
    else
    {
      result
    }
  }
}

stratified_resample<-function(log_weights)
{
  P = length(log_weights)
  W = exp(log_weights)
  cw = cumsum(W / sum(W))
  n = length(W)
  u = runif(n,0,1)/n
  v = u + (0:(n-1))/n
  
  j = 1
  indices = matrix(0,n)
  for (i in 1:n)
  {
    while(cw[j] < v[i])
    {
      j = j+1
    }
    indices[i] = j
  }
  
  indices
}

stratified_resample_u<-function(log_weights,u)
{
  P = length(log_weights)
  W = exp(log_weights)
  cw = cumsum(W / sum(W))
  n = length(W)
  u = u/n
  v = u + (0:(n-1))/n
  
  j = 1
  indices = matrix(0,n)
  for (i in 1:n)
  {
    while(cw[j] < v[i])
    {
      j = j+1
    }
    indices[i] = j
  }
  
  indices
}

# can take unnormalised weights
ess <- function(log_weights)
{
  the_ess = exp(2*logsumexp(log_weights) - logsumexp(2*log_weights))
  
  if (is.nan(the_ess))
  {
    0
  }
  else
  {
    the_ess
  }
}

particle_filter_initial_step = function(y,num_particles,inputs,initial_step_simulator,measurement_model_logevaluate,measurement_model_parameter_index,dim_x,resample_threshold=0.5)
{
  # Run initial step of PF.

  current_log_llhd_estimate = 0
  pf_states = list()[rep(1:num_particles)]
  pf_data = list()[rep(1:num_particles)]
  
  input_list = mclapply(as.list(c(1:num_particles)),function(x){c(inputs,x)})
  
  output = mclapply(input_list,function(th){initial_step_simulator(th)})
  pf_data = lapply(output, `[[`, 1)
  pf_states = lapply(output, `[[`, 2)
  
  pf_log_weights = matrix(-log(num_particles),num_particles)
  
  pf_log_weights = pf_log_weights + sapply(mclapply(pf_data,measurement_model_logevaluate,y,inputs[measurement_model_parameter_index]),c)
  
  log_sum_weights = logsumexp(pf_log_weights)
  current_log_llhd_estimate = current_log_llhd_estimate + log_sum_weights
  
  pf_log_weights = pf_log_weights-logsumexp(pf_log_weights)
  
  an_ess = ess(pf_log_weights)
  if (an_ess<resample_threshold*num_particles)
  {
    indices = stratified_resample(pf_log_weights)
    pf_states = pf_states[indices]
    pf_data = pf_data[indices]
    pf_log_weights = matrix(-log(num_particles),num_particles)
  }

  return(list(pf_data,pf_states,pf_log_weights,current_log_llhd_estimate))
}

particle_filter_step = function(y,num_particles,inputs,current_pf_states,current_pf_data,pf_log_weights,current_log_llhd_estimate,t,step_simulator,measurement_model_logevaluate,measurement_model_parameter_index,dim_x,resample_threshold=0.5)
{
  input_list = mclapply(as.list(c(1:num_particles)),function(x){c(inputs,x)})
  result = mcmapply(function(inp,sta){step_simulator(inp,sta,t)},input_list,current_pf_states)
  pf_data = lapply(1:num_particles,function(x){result[[1,x]]})
  pf_states = lapply(1:num_particles,function(x){result[[2,x]]})
  
  pf_log_weights = pf_log_weights + sapply(mclapply(pf_data,measurement_model_logevaluate,y,inputs[measurement_model_parameter_index]),c)
  
  log_sum_weights = logsumexp(pf_log_weights)
  current_log_llhd_estimate = current_log_llhd_estimate + log_sum_weights

  pf_log_weights = pf_log_weights-logsumexp(pf_log_weights)
  
  an_ess = ess(pf_log_weights)
  if (an_ess<resample_threshold*num_particles)
  {
    indices = stratified_resample(pf_log_weights)
    pf_states = pf_states[indices]
    pf_data = pf_data[indices]
    pf_log_weights = matrix(-log(num_particles),num_particles)
  }
  
  return(list(pf_data,pf_states,pf_log_weights,current_log_llhd_estimate))
}


particle_filter = function(y,num_particles,inputs,initial_step_simulator,step_simulator,measurement_model_logevaluate,measurement_model_parameter_index,dim_x,resample_threshold=0.5,termination_criterion=-Inf)
{
  # Run initial step of PF.
  
  if (is.matrix(y))
    T = nrow(y)
  else
    T = 1
  current_log_llhd_estimate = 0
  pf_states = list(list()[rep(1:num_particles)])[rep(1,T)]
  pf_data = list(list()[rep(1:num_particles)])[rep(1,T)]
  
  input_list = mclapply(as.list(c(1:num_particles)),function(x){c(inputs,x)})
  
  output = mclapply(input_list,function(th){initial_step_simulator(th)})
  pf_data[[1]] = lapply(output, `[[`, 1)
  pf_states[[1]] = lapply(output, `[[`, 2)
  
  pf_log_weights = matrix(-log(num_particles),num_particles)
  
  if (T>1)
    pf_log_weights = pf_log_weights + sapply(mclapply(pf_data[[1]],measurement_model_logevaluate,y[1,],inputs[measurement_model_parameter_index]),c)
  else
    pf_log_weights = pf_log_weights + sapply(mclapply(pf_data[[1]],measurement_model_logevaluate,y,inputs[measurement_model_parameter_index]),c)
  
  log_sum_weights = logsumexp(pf_log_weights)
  current_log_llhd_estimate = current_log_llhd_estimate + log_sum_weights
  
  if (current_log_llhd_estimate<termination_criterion)
  {
    return(list(current_log_llhd_estimate,matrix(0,0,0),matrix(0,0,0)))
  }
  
  pf_log_weights = pf_log_weights-logsumexp(pf_log_weights)
  
  an_ess = ess(pf_log_weights)
  if (an_ess<resample_threshold*num_particles)
  {
    indices = stratified_resample(pf_log_weights)
    pf_states[[1]] = pf_states[[1]][indices]
    pf_data[[1]] = pf_data[[1]][indices]
    pf_log_weights = matrix(-log(num_particles),num_particles)
  }
  
  # Run the rest of the PF.
  if (T>1)
  {
    for (t in 1:(T-1))
    {
      result = mcmapply(function(inp,sta){step_simulator(inp,sta,t+1)},input_list,pf_states[[t]])
      pf_data[[t+1]] = lapply(1:num_particles,function(x){result[[1,x]]})
      pf_states[[t+1]] = lapply(1:num_particles,function(x){result[[2,x]]})
      
      pf_log_weights = pf_log_weights + sapply(mclapply(pf_data[[t+1]],measurement_model_logevaluate,y[t+1,],inputs[measurement_model_parameter_index]),c)
      
      log_sum_weights = logsumexp(pf_log_weights)
      current_log_llhd_estimate = current_log_llhd_estimate + log_sum_weights
      if (current_log_llhd_estimate<termination_criterion)
      {
        return(list(current_log_llhd_estimate,matrix(0,0,0),matrix(0,0,0)))
      }
      pf_log_weights = pf_log_weights-logsumexp(pf_log_weights)
      
      an_ess = ess(pf_log_weights)
      if (an_ess<resample_threshold*num_particles)
      {
        indices = stratified_resample(pf_log_weights)
        pf_states = lapply(pf_states,`[`,indices)
        pf_data = lapply(pf_data,`[`,indices)
        pf_log_weights = matrix(-log(num_particles),num_particles)
      }
    }
  }
  
  return(list(current_log_llhd_estimate,pf_data,pf_states))
}


pmcmc_move = function(num_blocks,sample,current_inputs_for_simulator,index_of_input_in_param_order,index_of_param_in_input_order,proposals_index,proposals_simulate,num_within_block,priors_logevaluate,priors_index_in_proposal,priors_index,proposals_logevaluate,measurement_model_logevaluate,y,measurement_model_parameter_index,max_measurement_model,currentT,current_log_llhd_estimate,num_pmcmc_particles,initial_step_simulator,step_simulator,dim_x,pmcmc_resample_threshold,prop_covs)
{
  update_order = sample(1:num_blocks,num_blocks)
  for (a in 1:num_blocks)
  {
    b = update_order[a]
  
    proposed_parameter = sample
    proposed_inputs_for_simulator = current_inputs_for_simulator
    proposed_inputs_for_simulator[index_of_input_in_param_order] = proposed_parameter
    log_accept_prior = 0
    proposed_parameter[index_of_param_in_input_order[proposals_index[[b]]]] = proposals_simulate[[b]](conditioned_on_parameter=sample[index_of_param_in_input_order[proposals_index[[b]]]],other_parameters=sample,proposal_parameters=prop_covs[[b]])
    
    for (w in 1:num_within_block[b])
    {
      log_accept_prior = log_accept_prior + priors_logevaluate[[ priors_index_in_proposal[[b]][w] ]](proposed_parameter[ index_of_param_in_input_order[priors_index[[priors_index_in_proposal[[b]][w]]]] ]) - priors_logevaluate[[ priors_index_in_proposal[[b]][w] ]](sample[ index_of_param_in_input_order[priors_index[[priors_index_in_proposal[[b]][w]]]] ])
    }
    log_accept_prior = log_accept_prior + proposals_logevaluate[[b]](current_parameter=sample[index_of_param_in_input_order[proposals_index[[b]]]],conditioned_on_parameter=proposed_parameter[ index_of_param_in_input_order[proposals_index[[b]]] ],other_parameters=sample,proposal_parameters = prop_covs[[b]]) - proposals_logevaluate[[b]](current_parameter=proposed_parameter[ index_of_param_in_input_order[proposals_index[[b]]] ],conditioned_on_parameter=sample[ index_of_param_in_input_order[proposals_index[[b]]] ],other_parameters=sample,proposal_parameters = prop_covs[[b]])
    
    proposed_inputs_for_simulator[proposals_index[[b]]] = proposed_parameter[index_of_param_in_input_order[proposals_index[[b]]]]
      
    proposed_max_measurement_model = max_measurement_model(proposed_inputs_for_simulator[measurement_model_parameter_index])

    termination_criterion = log(runif(1)) - log_accept_prior + currentT*(-proposed_max_measurement_model) + current_log_llhd_estimate
    
    if (is.nan(termination_criterion))
      browser()
    
    result = particle_filter(y,num_pmcmc_particles,proposed_inputs_for_simulator,initial_step_simulator,step_simulator,measurement_model_logevaluate,measurement_model_parameter_index,dim_x,pmcmc_resample_threshold,termination_criterion)
    proposed_log_llhd_estimate = result[[1]]
    proposed_data = result[[2]]
    proposed_states = result[[3]]
    
    if (proposed_log_llhd_estimate>termination_criterion)
    {
      sample = proposed_parameter
      current_inputs_for_simulator = proposed_inputs_for_simulator
      current_log_llhd_estimate = proposed_log_llhd_estimate + currentT*max_measurement_model(current_inputs_for_simulator[measurement_model_parameter_index])
      states = proposed_states
      data = proposed_data
    }
    else
    {
      states = list()
      data = list()
    }
  }
  
  return(list(sample,data,states,current_log_llhd_estimate,current_inputs_for_simulator))
}

smc_squared = function(num_particles,num_pmcmc_particles,y,step_simulator,initial_step_simulator,measurement_model_logevaluate,max_measurement_model=function(x){0},inputs,priors_logevaluate,priors_simulate,priors_index=list(),proposals_simulate,proposals_logevaluate,proposals_index=list(),measurement_model_parameter_index,resample_threshold=0.5,pmcmc_resample_threshold=0.5,store_states=FALSE,store_data=TRUE,cores=1)
{
  options(mc.cores=cores)
  # Check dimensions of priors_evaluate,proposals_simulate,proposals_logevaluate,parameter_index
  
  # Case 1 - we have a single prior and single proposal, no inputs and no parameter index
  # In this case we need to set up parameter_index, enums_for_block, inputs, prior/proposal in block, and param count
  # case - we have single prior, and inputs and no param index
  if ( (length(priors_index)==0) || (length(proposals_index)==0) )
  {
    # No parameter_index supplied, so the assumption is that all of the dimensions in the priors, proposals, inputs and simulator refer to a vector in the same order.
    
    if ( (!is.function(priors_logevaluate)) || (!is.function(priors_simulate)) || (!is.function(proposals_simulate)) || (!is.function(proposals_logevaluate)) )
    {
      # throw error
      browser()
    }
    
    if ( (length(priors_index)>1) || (length(proposals_index)>1) ) {
      # throw error
      browser()
    }
    
    if (length(priors_index)==1) {
      if (length(inputs)!=length(priors_index[[1]])) {
        # throw error
        browser()
      }
      else {
        proposals_index = priors_index
      }
    } else if (length(proposals_index)==1) {
      if (length(inputs)!=length(proposals_index[[1]])) {
        # throw error
        browser()
      }
      else {
        priors_index = proposals_index
      }
    } else if (nrow(inputs)!=0) {
      priors_index = list()
      proposals_index = list()
      priors_index[[1]] = 1:length(inputs)
      proposals_index[[1]] = 1:length(inputs)
    } else {
      # throw error
      browser()
    }
  }
  
  # Sort out blocks...
  # If priors and proposals are not in blocks, then they need to be.
  
  if (!is.list(priors_index))
  {
    prior_input2 = priors_index
    priors_index = list()
    priors_index[[1]] = prior_input2
  }
  
  if (is.function(priors_logevaluate))
  {
    prior_input = priors_logevaluate
    priors_logevaluate = list()
    priors_logevaluate[[1]] = prior_input
    if ( is.list(priors_index) && (length(priors_index)!=1) )
    {
      # throw error
      browser()
    }
  }
  
  if (is.function(priors_simulate))
  {
    prior_input = priors_simulate
    priors_simulate = list()
    priors_simulate[[1]] = prior_input
    if ( is.list(priors_index) && (length(priors_index)!=1) )
    {
      # throw error
      browser()
    }
  }
  
  if (!is.list(proposals_index))
  {
    proposals_input2 = proposals_index
    proposals_index = list()
    proposals_index[[1]] = proposals_input2
  }
  
  if (is.function(proposals_logevaluate))
  {
    proposals_input = proposals_logevaluate
    proposals_logevaluate = list()
    proposals_logevaluate[[1]] = proposals_input
    if ( is.list(proposals_index) && (length(proposals_index)!=1) )
    {
      # throw error
      browser()
    }
  }
  
  if (is.function(proposals_simulate))
  {
    proposals_input = proposals_simulate
    proposals_simulate = list()
    proposals_simulate[[1]] = proposals_input
    if ( is.list(proposals_index) && (length(proposals_index)!=1) )
    {
      # throw error
      browser()
    }
  }
  
  # Also need to check indices are lists!
  
  if (length(priors_index)!=length(priors_logevaluate))
  {
    # throw error
    browser()
  }
  
  if (length(proposals_index)!=length(proposals_logevaluate))
  {
    # throw error
    browser()
  }
  
  if (length(proposals_index)!=length(proposals_simulate))
  {
    # throw error
    browser()
  }
  
  # Give num_blocks.
  num_blocks = length(proposals_logevaluate)
  
  # Check that no block splits a prior.
  for (i in 1:length(priors_index))
  {
    split = FALSE
    for (j in 1:length(proposals_index))
    {
      inter = intersect(priors_index[[i]],proposals_index[[j]])
      if (length(inter)<length(priors_index[[i]]))
      {
        if (length(inter)!=0)
          split = TRUE
      }
    }
    if (split==TRUE)
    {
      # throw error
      browser()
    }
  }
  
  num_within_block = matrix(0,num_blocks)
  
  # Each block gets a list of priors to use.
  priors_index_in_proposal = list()
  
  for (j in 1:length(proposals_index))
  {
    current_priors_index = c()
    for (i in 1:length(priors_index))
    {
      inter = intersect(priors_index[[i]],proposals_index[[j]])
      if (length(inter)!=0)
      {
        current_priors_index = c(current_priors_index,i)
      }
    }
    priors_index_in_proposal[[j]] = current_priors_index
    num_within_block[j] = length(current_priors_index)
  }
  
  max_index = 0
  num_params = 0
  index_of_input_in_param_order = c()
  index_of_param_in_input_order = rep(0,length(inputs))
  for (i in 1:length(priors_index))
  {
    index_of_param_in_input_order[priors_index[[i]]] = (num_params+1):(num_params+length(priors_index[[i]]))
    num_params = num_params + length(priors_index[[i]])
    index_of_input_in_param_order = c(index_of_input_in_param_order,priors_index[[i]])
    
    for (k in 1:length(priors_index[[i]]))
    {
      if (priors_index[[i]][k]>max_index)
      {
        max_index = priors_index[[i]][k]
      }
    }
  }
  
  if (nrow(inputs)==0) {
    inputs = matrix(0,max_index)
  }
  else {
    if (max_index>length(inputs))
    {
      # throw error
      browser()
    }
  }
  
  num_params = length(index_of_input_in_param_order)
  
  current_max_measurement_model = max_measurement_model(inputs[measurement_model_parameter_index])
  
  # Check simulator doesn't throw error.
  # Check dimension of simulator output matches the real data.
  sim_data = matrix(0,0,0)
  sim_states = matrix(0,0,0)
  example_output = initial_step_simulator(c(inputs,0))
  
  if (is.list(example_output))
  {
    state_output = TRUE
  }
  else
  {
    state_output = FALSE
  }
  
  if (state_output==TRUE)
  {
    if (store_states==TRUE)
    {
      sim_states = example_output[[2]]
    }
    
    sim_data = example_output[[1]]
  }
  else
  {
    # throw error
    browser()
  }
  
  # Check that both the data and the sim_data work in the measurement model.
  ##########################################################################
  # Not yet implemented!
  
  dim_x = length(sim_data)
  T = nrow(y)
  
  print(1)
  
  particles = matrix(0,num_particles,num_params)
  
  current_log_llhd_estimates = matrix(0,num_particles)
  
  proposed_inputs_for_simulator = t(matrix(rep(inputs,num_particles),length(inputs),num_particles))
  
  for (j in 1:length(priors_index))
  {
    particles[,index_of_param_in_input_order[priors_index[[j]]]] = priors_simulate[[j]](num_particles)
    proposed_inputs_for_simulator[,priors_index[[j]]] = particles[,index_of_param_in_input_order[priors_index[[j]]]]
  }
  
  current_inputs_for_simulator = proposed_inputs_for_simulator
  
  if (is.matrix(proposed_inputs_for_simulator))
    rows = lapply(1:num_particles,function(i){proposed_inputs_for_simulator[i,]})
  else
    rows = lapply(1:num_particles,function(i){proposed_inputs_for_simulator[i]})
  
  pf_states = list(list(list()[rep(1:num_pmcmc_particles)])[rep(1,num_particles)])[rep(1,T)]
  pf_data = list(list(list()[rep(1:num_pmcmc_particles)])[rep(1,num_particles)])[rep(1,T)]
  pf_log_weights = matrix(0,num_pmcmc_particles,num_particles)
  pf_weights = list(list()[rep(1,num_particles)])[rep(1,T)]
  
  current_w = matrix(-log(num_particles),num_particles)
  
  pf_output = lapply(rows,function(th){particle_filter_initial_step(y[1,],num_pmcmc_particles,th,initial_step_simulator,measurement_model_logevaluate,measurement_model_parameter_index,dim_x,pmcmc_resample_threshold)})
  
  pf_data[[1]] = lapply(pf_output, `[[`, 1)
  pf_states[[1]] = lapply(pf_output, `[[`, 2)
  pf_log_weights = matrix(unlist(lapply(pf_output, `[[`, 3)),num_pmcmc_particles,num_particles)
  current_log_llhd_estimates = unlist(lapply(pf_output, `[[`, 4)) + apply(as.matrix(current_inputs_for_simulator[,measurement_model_parameter_index]),1,max_measurement_model)
  
  current_w = current_w + current_log_llhd_estimates
  
  an_ess = ess(current_w)
  print(an_ess)
  if (an_ess<resample_threshold*num_particles)
  {
    indices = stratified_resample(current_w)
    pf_states[[1]] = pf_states[[1]][indices]
    pf_data[[1]] = pf_data[[1]][indices]
    pf_log_weights = pf_log_weights[,indices]
    current_w = matrix(-log(num_particles),num_particles)
    particles = particles[indices,]
    current_log_llhd_estimates = current_log_llhd_estimates[indices]
    
    if (is.matrix(particles))
    {
      particle_list = lapply(1:num_particles,function(i){particles[i,]})
    }
    else
    {
      particle_list = lapply(1:num_particles,function(i){particles[i]})
    }
    if (is.matrix(current_inputs_for_simulator))
      input_list = lapply(1:num_particles,function(i){current_inputs_for_simulator[i,]})
    else
      input_list = lapply(1:num_particles,function(i){current_inputs_for_simulator[i]})
    log_llhd_list = lapply(1:num_particles,function(i){current_log_llhd_estimates[i]})
    
    prop_covs = list()
      
    for (b in 1:num_blocks)
    {
      if (length(proposals_index[[b]])==num_params)
      {
        if (is.matrix(particles))
          prop_covs[[b]] = make.positive.definite(cov(particles))
        else
          prop_covs[[b]] = var(particles)
      }
      else
      {
        current_indices = index_of_param_in_input_order[proposals_index[[b]]]
        conditioned_on_indices = setdiff(1:num_params,current_indices)
        Z = as.data.frame(particles[,conditioned_on_indices])
        fit = lm(particles[,current_indices]~.,data=Z)
        if (is.matrix(particles[,current_indices]))
          prop_covs[[b]] = make.positive.definite(cov(residuals(fit)))
        else
          prop_covs[[b]] = var(residuals(fit))
      }
    }
    
    pf_output = mlapply(list(particle_list,input_list,log_llhd_list),function(a,b,e){pmcmc_move(num_blocks,a,b,index_of_input_in_param_order,index_of_param_in_input_order,proposals_index,proposals_simulate,num_within_block,priors_logevaluate,priors_index_in_proposal,priors_index,proposals_logevaluate,measurement_model_logevaluate,y[1,],measurement_model_parameter_index,max_measurement_model,1,e,num_pmcmc_particles,initial_step_simulator,step_simulator,dim_x,pmcmc_resample_threshold,prop_covs)})

    particles = t(matrix(unlist(lapply(pf_output, `[[`, 1)),num_params,num_particles))

    pf_data_out = lapply(pf_output, `[[`, 2)
    which_changed = which(apply(as.matrix(1:num_particles),1,function(i){length(unlist(pf_data_out[[i]]))>0}))
    print(length(which_changed))
    if (length(which_changed)>0)
    {
      pf_data[[1]][which_changed] = lapply(1:length(which_changed),function(i){pf_data_out[[which_changed[i]]][[1]]})
    }

    pf_states_out = lapply(pf_output, `[[`, 3)
    which_changed = which(apply(as.matrix(1:num_particles),1,function(i){length(unlist(pf_states_out[[i]]))>0}))
    if (length(which_changed)>0)
    {
      pf_states[[1]][which_changed] = lapply(1:length(which_changed),function(i){pf_states_out[[which_changed[i]]][[1]]})
    }
    current_log_llhd_estimates = unlist(lapply(pf_output, `[[`, 4))
    current_inputs = t(matrix(unlist(lapply(pf_output, `[[`, 5)),length(inputs),num_particles))
    
  }
  
  # Other steps.
  
  for (t in 2:T)
  {
    
    print(t)
    
    weights_list = lapply(1:num_particles,function(i){pf_log_weights[,i]})
    log_llhd_list = lapply(1:num_particles,function(i){current_log_llhd_estimates[i]})
    input_list = lapply(1:num_particles,function(i){current_inputs_for_simulator[i,]})
    
    pf_output = mlapply(list(input_list,pf_states[[t-1]],pf_data[[t-1]],weights_list,log_llhd_list),function(a,b,c,d,e){particle_filter_step(y[t,],num_pmcmc_particles,a,b,c,d,e,t,step_simulator,measurement_model_logevaluate,measurement_model_parameter_index,dim_x,pmcmc_resample_threshold)})
    
    previous_log_llhd_estimates = current_log_llhd_estimates
    
    pf_data[[t]] = lapply(pf_output, `[[`, 1)
    pf_states[[t]] = lapply(pf_output, `[[`, 2)
    
    pf_log_weights = matrix(unlist(lapply(pf_output, `[[`, 3)),num_pmcmc_particles,num_particles)
    current_log_llhd_estimates = unlist(lapply(pf_output, `[[`, 4)) + apply(as.matrix(current_inputs_for_simulator[,measurement_model_parameter_index]),1,max_measurement_model)
    
    current_w = current_w + current_log_llhd_estimates - previous_log_llhd_estimates
    
    an_ess = ess(current_w)
    print(an_ess)
    if (an_ess<resample_threshold*num_particles)
    {
      indices = stratified_resample(current_w)
      for (k in 1:t)
      {
        pf_states[[k]] = pf_states[[k]][indices]
        pf_data[[k]] = pf_data[[k]][indices]
      }
      
      pf_log_weights = pf_log_weights[,indices]
      current_w = matrix(-log(num_particles),num_particles)
      particles = particles[indices,]
      current_log_llhd_estimates = current_log_llhd_estimates[indices]
      current_inputs_for_simulator = current_inputs_for_simulator[indices,]
      
      if (is.matrix(particles))
      {
        particle_list = lapply(1:num_particles,function(i){particles[i,]})
      }
      else
      {
        particle_list = lapply(1:num_particles,function(i){particles[i]})
      }
      if (is.matrix(current_inputs_for_simulator))
        input_list = lapply(1:num_particles,function(i){current_inputs_for_simulator[i,]})
      else
        input_list = lapply(1:num_particles,function(i){current_inputs_for_simulator[i]})
      log_llhd_list = lapply(1:num_particles,function(i){current_log_llhd_estimates[i]})
      
      prop_covs = list()
      
      for (b in 1:num_blocks)
      {
        if (length(proposals_index[[b]])==num_params)
        {
          if (is.matrix(particles))
            prop_covs[[b]] = make.positive.definite(cov(particles))
          else
            prop_covs[[b]] = var(particles)
        }
        else
        {
          current_indices = index_of_param_in_input_order[proposals_index[[b]]]
          conditioned_on_indices = setdiff(1:num_params,current_indices)
          Z = as.data.frame(particles[,conditioned_on_indices])
          fit = lm(particles[,current_indices]~.,data=Z)
          if (is.matrix(particles[,current_indices]))
            prop_covs[[b]] = make.positive.definite(cov(residuals(fit)))
          else
            prop_covs[[b]] = var(residuals(fit))
        }
      }
      
      pf_output = mlapply(list(particle_list,input_list,log_llhd_list),function(a,b,e){pmcmc_move(num_blocks,a,b,index_of_input_in_param_order,index_of_param_in_input_order,proposals_index,proposals_simulate,num_within_block,priors_logevaluate,priors_index_in_proposal,priors_index,proposals_logevaluate,measurement_model_logevaluate,y[1:t,],measurement_model_parameter_index,max_measurement_model,t,e,num_pmcmc_particles,initial_step_simulator,step_simulator,dim_x,pmcmc_resample_threshold,prop_covs)})
      
      particles = t(matrix(unlist(lapply(pf_output, `[[`, 1)),num_params,num_particles))
      pf_data_out = lapply(pf_output, `[[`, 2)
      which_changed = which(apply(as.matrix(1:num_particles),1,function(i){length(unlist(pf_data_out[[i]]))>0}))
      print(length(which_changed))
      if (length(which_changed)>0)
      {
        for (k in 1:t)
          pf_data[[k]][which_changed] = lapply(1:length(which_changed),function(i){pf_data_out[[which_changed[i]]][[k]]})
      }
      pf_states_out = lapply(pf_output, `[[`, 3)
      which_changed = which(apply(as.matrix(1:num_particles),1,function(i){length(unlist(pf_states_out[[i]]))>0}))
      if (length(which_changed)>0)
      {
        for (k in 1:t)
          pf_states[[k]][which_changed] = lapply(1:length(which_changed),function(i){pf_states_out[[which_changed[i]]][[k]]})
      }
      
      current_log_llhd_estimates = unlist(lapply(pf_output, `[[`, 4))
      
    }
  }

  
  if ((store_data=TRUE) && (store_states==TRUE))
  {
    return(list(particles,current_w,index_of_input_in_param_order,pf_data,pf_states,pf_log_weights))
  } else if ((store_data=TRUE) && (store_states==FALSE))
  {
    return(list(particles,current_w,index_of_input_in_param_order,pf_data))
  } else if ((store_data=FALSE) && (store_states==TRUE))
  {
    return(list(particles,current_w,index_of_input_in_param_order,pf_states))
  } else if ((store_data=FALSE) && (store_states==FALSE))
  {
    return(list(particles,current_w,index_of_input_in_param_order))
  }
}