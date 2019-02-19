library(LaplacesDemon)
library(gsubfn)
require(parallel)

is.string = function(x)
{
  is.character(x) && length(x) == 1
}

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

particle_filter = function(y,num_particles,inputs,initial_step_simulator,step_simulator,measurement_model_logevaluate,measurement_model_parameter_index,dim_x,resample_threshold=0.5,termination_criterion=-Inf)
{
  # Run initial step of PF.
  
  if (is.matrix(y))
    T = nrow(y)
  else
    T = 1
  current_log_llhd_estimate = 0
  pf_states = list(list()[rep(1:num_particles)])[rep(1,T)]
  pf_data = list(list()[rep(1:num_particles)])[rep(1,T)]#array(0,dim=c(dim_x,num_particles,T))
  
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
  
  if (is.infinite(current_log_llhd_estimate))
  {
    return(list(current_log_llhd_estimate,matrix(0,0,0),matrix(0,0,0)))
  }

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
  for (t in 1:(T-1))
  {
    result = mcmapply(function(inp,sta){step_simulator(inp,sta,t+1)},input_list,pf_states[[t]])
    pf_data[[t+1]] = lapply(1:num_particles,function(x){result[[1,x]]})
    pf_states[[t+1]] = lapply(1:num_particles,function(x){result[[2,x]]})
    
    pf_log_weights = pf_log_weights + sapply(mclapply(pf_data[[t+1]],measurement_model_logevaluate,y[t+1,],inputs[measurement_model_parameter_index]),c)
    
    log_sum_weights = logsumexp(pf_log_weights)
    current_log_llhd_estimate = current_log_llhd_estimate + log_sum_weights
    if (is.infinite(current_log_llhd_estimate))
    {
      return(list(current_log_llhd_estimate,matrix(0,0,0),matrix(0,0,0)))
    }

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
  
  # Randomly sample one index to give the state to keep.
  index = sample(x = 1:num_particles, 1, replace = TRUE, prob = exp(pf_log_weights))
  current_state = lapply(pf_states,`[`,index)
  current_data = lapply(pf_data,`[`,index)
  current_data = t(matrix(unlist(lapply(current_data, `[[`, 1)),dim_x,T))
  
  return(list(current_log_llhd_estimate,current_data,current_state))
}

particle_mcmc = function(num_mcmc,y,step_simulator,initial_step_simulator,measurement_model_logevaluate,max_measurement_model=function(x){0},inputs,priors_logevaluate,priors_index=list(),proposals_simulate,proposals_logevaluate,proposals_index=list(),measurement_model_parameter_index,num_particles,resample_threshold=0.5,store_states=FALSE,store_data=TRUE,cores=1)
{
  options(mc.cores=cores)
  # Check dimensions of priors_evaluate,proposals_simulate,proposals_logevaluate,parameter_index
  
  # Case 1 - we have a single prior and single proposal, no inputs and no parameter index
  # In this case we need to set up parameter_index, enums_for_block, inputs, prior/proposal in block, and param count
  # case - we have single prior, and inputs and no param index
  if ( (length(priors_index)==0) || (length(proposals_index)==0) )
  {
    # No parameter_index supplied, so the assumption is that all of the dimensions in the priors, proposals, inputs and simulator refer to a vector in the same order.
    
    if ( (!is.function(priors_logevaluate)) || (!is.function(proposals_simulate)) || (!is.function(proposals_logevaluate)) )
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
    } else if (length(inputs)!=0) {
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
  
  #browser()
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
  
  if (length(inputs)==0) {
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
  
  logu = matrix(log(runif(num_mcmc*num_blocks)),num_blocks,num_mcmc)
  
  sample = matrix(0,num_params,num_mcmc)
  
  # Run initial step of PF.
  result = particle_filter(y,num_particles,inputs,initial_step_simulator,step_simulator,measurement_model_logevaluate,measurement_model_parameter_index,dim_x,resample_threshold)
  current_log_llhd_estimate = result[[1]]
  current_data = result[[2]]
  current_state = result[[3]]
  
  current_log_llhd_estimate = current_log_llhd_estimate + T*max_measurement_model(inputs[measurement_model_parameter_index])

  if (store_states==TRUE)
  {
    states = list()[rep(1,num_mcmc)]
    states[[1]] = current_state
  }
  
  if (store_data==TRUE)
  {
    data = array(0,dim=c(T,dim_x,num_mcmc))
    if (length(current_data)>0)
      data[,,1] = current_data
    else
      data[,,1] = matrix(0,T,dim_x)
  }
  else
    data = matrix(0,0,0)
  
  # if (is.infinite(current_log_llhd_estimate))
  # {
  #   # throw error
  #   browser()
  # }
  
  sample[,1] = inputs[index_of_input_in_param_order]
  
  current_inputs_for_simulator = inputs
  
  for (i in 2:num_mcmc)
  {
    if (i%%1==0)
    {
      print(sprintf("Current iteration: %i", i))
    }
    sample[,i] = sample[,i-1]
    if (store_states==TRUE)
      states[[i]] = states[[i-1]]
    if (store_data==TRUE)
      data[,,i] = data[,,i-1]
    
    update_order = sample(1:num_blocks,num_blocks)
    for (a in 1:num_blocks)
    {
      b = update_order[a]
      
      proposed_parameter = sample[,i]
      proposed_inputs_for_simulator = current_inputs_for_simulator
      proposed_inputs_for_simulator[index_of_input_in_param_order] = proposed_parameter
      log_accept_prior = 0
      proposed_parameter[index_of_param_in_input_order[proposals_index[[b]]]] = proposals_simulate[[b]](conditioned_on_parameter=sample[index_of_param_in_input_order[proposals_index[[b]]],i],other_parameters=sample[,i])
      
      for (w in 1:num_within_block[b])
      {
        log_accept_prior = log_accept_prior + priors_logevaluate[[ priors_index_in_proposal[[b]][w] ]](proposed_parameter[ index_of_param_in_input_order[priors_index[[priors_index_in_proposal[[b]][w]]]] ]) - priors_logevaluate[[ priors_index_in_proposal[[b]][w] ]](sample[ index_of_param_in_input_order[priors_index[[priors_index_in_proposal[[b]][w]]]] ,i])
      }
      log_accept_prior = log_accept_prior + proposals_logevaluate[[b]](current_parameter=sample[index_of_param_in_input_order[proposals_index[[b]]],i],conditioned_on_parameter=proposed_parameter[ index_of_param_in_input_order[proposals_index[[b]]] ],other_parameters=sample[,i]) - proposals_logevaluate[[b]](current_parameter=proposed_parameter[ index_of_param_in_input_order[proposals_index[[b]]] ],conditioned_on_parameter=sample[ index_of_param_in_input_order[proposals_index[[b]]] ,i],other_parameters=sample[,i])
      
      proposed_inputs_for_simulator[proposals_index[[b]]] = proposed_parameter[index_of_param_in_input_order[proposals_index[[b]]]]
      
      proposed_max_measurement_model = max_measurement_model(proposed_inputs_for_simulator[measurement_model_parameter_index])
      
      termination_criterion = logu[b,i] - log_accept_prior + T*(-proposed_max_measurement_model) + current_log_llhd_estimate
      
      result = particle_filter(y,num_particles,proposed_inputs_for_simulator,initial_step_simulator,step_simulator,measurement_model_logevaluate,measurement_model_parameter_index,dim_x,resample_threshold,termination_criterion)
      proposed_log_llhd_estimate = result[[1]]
      proposed_data = result[[2]]
      proposed_states = result[[3]]
      
      if (!is.infinite(proposed_log_llhd_estimate))
      {
        if (proposed_log_llhd_estimate>termination_criterion)
        {
          sample[,i] = proposed_parameter
          current_inputs_for_simulator = proposed_inputs_for_simulator
          current_log_llhd_estimate = proposed_log_llhd_estimate + T*max_measurement_model(current_inputs_for_simulator[measurement_model_parameter_index])
          print("Accepted")
          if (store_states==TRUE)
            states[[i]] = proposed_states
          if (store_data==TRUE)
            data[,,i] = proposed_data
        }
      }
    }
  }
  
  if ((store_data=TRUE) && (store_states==TRUE))
  {
    return(list(sample,index_of_input_in_param_order,data,states))
  } else if ((store_data=TRUE) && (store_states==FALSE))
  {
    return(list(sample,index_of_input_in_param_order,data))
  } else if ((store_data=FALSE) && (store_states==TRUE))
  {
    return(list(sample,index_of_input_in_param_order,states))
  } else if ((store_data=FALSE) && (store_states==FALSE))
  {
    return(list(sample,index_of_input_in_param_order))
  }
}