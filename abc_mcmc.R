library(LaplacesDemon)
#library(rapportools)
library(gsubfn)
#library(svMisc)

# Want to add:
# - adaptive MCMC that uses ESS
# - r-hit
# - non-uniform kernel
# - more than one ABC sim

is.string = function(x)
{
  is.character(x) && length(x) == 1
}

abc_mcmc = function(num_mcmc,y,simulator,inputs,priors_logevaluate,priors_index=list(),proposals_simulate,proposals_logevaluate,proposals_index=list(),tolerance=Inf,summary_scaling=1,num_abc_sims=1,abc_summary="identity",abc_kernel="uniform",distance="euclidean",store_states=FALSE,store_data=TRUE,time_initialising=2)
{
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
  
  # Check simulator doesn't throw error.
  # Check dimension of simulator output matches the real data.
  sim_data = matrix(0,0,0)
  sim_states = matrix(0,0,0)
  example_output = simulator(inputs)
  
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
      num_states = length(states)
    }
    
    sim_data = example_output[[1]]
  }
  else
  {
    sim_data = example_output
  }
  
  if (length(sim_data)!=length(y))
  {
    # throw error
    browser()
  }
  
  num_data = length(sim_data)
  
  if (is.string(abc_summary))
  {
    if (abc_summary=="identity") {
      abc_summary = function(data){return(data)}
    } else {
      # throw error
      browser()
    }
  } else {
    if (!is.function(abc_summary))
    {
      # throw error
      browser()
    }
  }
  
  if (is.string(abc_kernel))
  {
    if (!(abc_kernel=="uniform")) {
      # throw error
      browser()
    }
  } else {
    if (!is.function(abc_kernel))
    {
      # throw error
      browser()
    }
  }
  
  if (is.string(distance))
  {
    if (distance=="euclidean") {
      distance = function(observed,simulated) {
        return(sqrt(sum(sum(((simulated-observed)*summary_scaling)^2))))
      }
    } else {
      # throw error
      browser()
    }
  } else {
    if (!is.function(distance))
    {
      # throw error
      browser()
    }
  }
  
  sum_y = abc_summary(y)
  
  # Must be a simulation for initial conditions that passes the tolerance.
  # Try for 2 minutes default.
  
  if (tolerance<Inf)
  {
    sum_x = abc_summary(sim_data)
    start_clock_time <- Sys.time()
    min_dist = distance(sum_x,sum_y)
    while (distance(sum_x,sum_y)>tolerance)
    {
      sim_data = simulator(inputs)
      
      if (state_output==TRUE)
      {
        if (store_states==TRUE)
        {
          sim_states = sim_data[[2]]
        }
        
        sim_data = sim_data[[1]]
      }
      
      sum_x = abc_summary(sim_data)
      min_dist = min(c(distance(sum_x,sum_y),min_dist))
      
      current_clock_time <- Sys.time()
      difference = unclass(current_clock_time-start_clock_time)
      if ( (attributes(difference)$units=="mins") && (difference>time_initialising) )
      {
        # throw error
        print(min_dist)
        browser()
      }
    }
  } else {
    # Simulate for some time.
    # Take epsilon to be the minimum distance.
  }
  
  logu = matrix(log(runif(num_mcmc*num_blocks)),num_blocks,num_mcmc)
  
  sample = matrix(0,num_params,num_mcmc)
  if (store_states==TRUE)
  {
    states = matrix(0,num_states,num_mcmc)  
    states[,1] = matrix(sim_states,length(sim_states))
  }
  else
    states = matrix(0,0,0)
  if (store_data==TRUE)
  {
    data = matrix(0,num_data,num_mcmc)
    data[,1] = matrix(sim_data,length(sim_data))
  }
  else
    data = matrix(0,0,0)
  
  sample[,1] = inputs[index_of_input_in_param_order]
  
  current_inputs_for_simulator = inputs
  
  for (i in 2:num_mcmc)
  {
    if (i%%1000==0)
      print(sprintf("Current iteration: %i", i))
    sample[,i] = sample[,i-1]
    if (store_states==TRUE)
      states[,i] = states[,i-1]
    if (store_data==TRUE)
      data[,i] = data[,i-1]
    
    update_order = sample(1:num_blocks,num_blocks)
    for (a in 1:num_blocks)
    {
      b = update_order[a]
      
      stored_dist = c()
      
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
      
      if (logu[b,i]<log_accept_prior)
      {
        proposed_inputs_for_simulator[proposals_index[[b]]] = proposed_parameter[index_of_param_in_input_order[proposals_index[[b]]]]
        proposed_data = simulator(proposed_inputs_for_simulator)
        if (state_output==TRUE)
        {
          if (store_states==TRUE)
          {
            proposed_states = proposed_data[[2]]
          }
          
          proposed_data = proposed_data[[1]]
        }
        
        sum_x = abc_summary(proposed_data)
        if (distance(sum_x,sum_y)<=tolerance)
        {
          sample[,i] = proposed_parameter
          current_inputs_for_simulator = proposed_inputs_for_simulator
          if (store_states==TRUE)
            states[,i] = proposed_states
          if (store_data==TRUE)
            data[,i] = proposed_data
        }
      }
      stored_dist = c(stored_dist,distance(sum_x,sum_y))
    }
  }
  
  return(list(sample,data,states,index_of_input_in_param_order))
}