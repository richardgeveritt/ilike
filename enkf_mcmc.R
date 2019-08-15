library(LaplacesDemon)
library(gsubfn)
#library(svMisc)
library(mvtnorm)

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

enkf = function(y,num_particles,inputs,initial_step_simulator,step_simulator,H,measurement_model_parameter_index,dim_x,resample_threshold=0.5,termination_criterion=-Inf,shift=TRUE)
{
  # Run initial step of EnKF.
  meas_dim = nrow(H)
  if (length(measurement_model_parameter_index)==1)
    R = diag(inputs[measurement_model_parameter_index]^2,meas_dim)
  else
    R = diag(inputs[measurement_model_parameter_index]^2)
  max_measurement_model = dmvnorm(t(matrix(0,meas_dim)), t(matrix(0,meas_dim)), R, log=TRUE)

  if (is.matrix(y))
    T = nrow(y)
  else
    T = 1
  current_log_llhd_estimate = 0
  pf_states = list(list()[rep(1:num_particles)])[rep(1,T)]
  pf_data = list(list()[rep(1:num_particles)])[rep(1,T)]#array(0,dim=c(dim_x,num_particles,T))
  
  input_list = lapply(as.list(c(1:num_particles)),function(x){c(inputs,x)})
  
  output = lapply(input_list,function(th){initial_step_simulator(th)})
  pf_data[[1]] = lapply(output, `[[`, 1)
  pf_states[[1]] = lapply(output, `[[`, 2)

  unlisted_state = unlist(pf_states[[1]])
  raw_state_vec = matrix(unlisted_state,length(pf_states[[1]][[1]]),length(unlisted_state)/length(pf_states[[1]][[1]]))
  state_vec = matrix(0,nrow(raw_state_vec),ncol(raw_state_vec))
  
  keptx = unique(which((!is.na(raw_state_vec))&(!is.infinite(raw_state_vec)),arr.ind = TRUE)[,1])
  kepty = unique(which((!is.na(raw_state_vec))&(!is.infinite(raw_state_vec)),arr.ind = TRUE)[,2])
  
  state_vec[keptx,kepty] = raw_state_vec[keptx,kepty]

  mu_tilde = rowMeans(state_vec)
  Sigma_tilde = cov(t(state_vec))
  
  latter_part = Sigma_tilde%*%t(H)
  latter_part[is.nan(latter_part)] = 0
  HSH = H%*%latter_part
  
  Hbymu = H%*%mu_tilde
  Hbymu[is.nan(Hbymu)] = 0
  
  if (T>1)
    current_log_llhd_estimate = current_log_llhd_estimate + dmvnorm(y[1,], Hbymu, HSH + R,log=TRUE) - max_measurement_model
  else
    current_log_llhd_estimate = current_log_llhd_estimate + dmvnorm(y, Hbymu, HSH + R,log=TRUE) - max_measurement_model
  
  if (is.infinite(current_log_llhd_estimate))
  {
    return(list(current_log_llhd_estimate,matrix(0,0,0),matrix(0,0,0)))
  }

  if (current_log_llhd_estimate<termination_criterion)
  {
    return(list(current_log_llhd_estimate,matrix(0,0,0),matrix(0,0,0)))
  }
  
  K = latter_part%*%solve(HSH + R)
  
  state_store = list()
  data_store = matrix(0,nrow(y),ncol(y))
  
  if (shift==FALSE)
  {
    if (T>1)
      mu_hat = mu_tilde + K%*%(y[1,] - Hbymu)
    else
      mu_hat = mu_tilde + K%*%(y - Hbymu)
    
    Sigma_hat = (diag(1,nrow(raw_state_vec)) - K%*%H)%*%Sigma_tilde
    Sigma_hat = as.matrix(forceSymmetric(Sigma_hat))
    
    state_store[[1]]= mu_hat
    
    data_store[1,] = H%*%mu_hat
  }
  
  # Run the rest of the EnKF.
  for (t in 1:(T-1))
  {
    if (shift==TRUE)
    {
      pf_states[[t]] = lapply(pf_states[[t]],function(x){state = x; state[is.infinite(state)] = 0; Hx = H%*%state; y_p = rmvnorm(1,Hx,R); x+K%*%(y[t,]-y_p)})
    }
    else {
      pf_states[[t]] = split(rmvnorm(num_particles,mu_hat,Sigma_hat),row(t(state_vec)))
    }
    
    result = mapply(function(inp,sta){step_simulator(inp,sta,t+1)},input_list,pf_states[[t]])
    pf_data[[t+1]] = lapply(1:num_particles,function(x){result[[1,x]]})
    pf_states[[t+1]] = lapply(1:num_particles,function(x){result[[2,x]]})
    
    unlisted_state = unlist(pf_states[[t+1]])  
    raw_state_vec = matrix(unlisted_state,length(pf_states[[t+1]][[1]]),length(unlisted_state)/length(pf_states[[t+1]][[1]]))
    state_vec = matrix(0,nrow(raw_state_vec),ncol(raw_state_vec))
    
    keptx = unique(which((!is.na(raw_state_vec))&(!is.infinite(raw_state_vec)),arr.ind = TRUE)[,1])
    kepty = unique(which((!is.na(raw_state_vec))&(!is.infinite(raw_state_vec)),arr.ind = TRUE)[,2])
    
    state_vec[keptx,kepty] = raw_state_vec[keptx,kepty]
    
    mu_tilde = rowMeans(state_vec)
    Sigma_tilde = cov(t(state_vec))
    if (any(is.nan(Sigma_tilde)))
      browser()
    
    latter_part = Sigma_tilde%*%t(H)
    latter_part[is.nan(latter_part)] = 0
    
    HSH = H%*%latter_part
    
    Hbymu = H%*%mu_tilde
    Hbymu[is.nan(Hbymu)] = 0
    
    if(is.na(dmvnorm(y[t+1,], Hbymu, HSH + R,log=TRUE)))
      browser()
    current_log_llhd_estimate = current_log_llhd_estimate + dmvnorm(y[t+1,], Hbymu, HSH + R,log=TRUE) - max_measurement_model
    
    if (is.infinite(current_log_llhd_estimate))
    {
      return(list(current_log_llhd_estimate,matrix(0,0,0),matrix(0,0,0)))
    }

    if (current_log_llhd_estimate<termination_criterion)
    {
      return(list(current_log_llhd_estimate,matrix(0,0,0),matrix(0,0,0)))
    }
    
    K = latter_part%*%solve(HSH + R)
    
    if (shift==FALSE)
    {
      mu_hat = mu_tilde + K%*%(y[t+1,] - H%*%mu_tilde)
      Sigma_hat = (diag(1,nrow(raw_state_vec)) - K%*%H)%*%Sigma_tilde
      Sigma_hat = as.matrix(forceSymmetric(Sigma_hat))
      if (any(is.nan(Sigma_hat)))
        return(list(-Inf,matrix(0,0,0),matrix(0,0,0)))
      
      state_store[[t+1]]= mu_hat
      data_store[t+1,] = H%*%mu_hat
    }
  }
  
  return(list(current_log_llhd_estimate,data_store,state_store))
}

enkf_mcmc = function(num_mcmc,y,step_simulator,initial_step_simulator,H,inputs,priors_logevaluate,priors_index=list(),proposals_simulate,proposals_logevaluate,proposals_index=list(),measurement_model_parameter_index,num_particles,store_states=FALSE,store_data=TRUE,cores=1)
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
  
  meas_dim = nrow(H)
  if (length(measurement_model_parameter_index)==1)
    R = diag(inputs[measurement_model_parameter_index]^2,meas_dim)
  else
    R = diag(inputs[measurement_model_parameter_index]^2)

  current_max_measurement_model = dmvnorm(t(matrix(0,meas_dim)), t(matrix(0,meas_dim)), R, log=TRUE)

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
  
  # Run initial step of EnKF.
  
  result = enkf(y,num_particles,inputs,initial_step_simulator,step_simulator,H,measurement_model_parameter_index,dim_x,resample_threshold)
  current_log_llhd_estimate = result[[1]]
  current_data = result[[2]]
  current_state = result[[3]]
  
  if (length(measurement_model_parameter_index)==1)
    R = diag(inputs[measurement_model_parameter_index]^2,meas_dim)
  else
    R = diag(inputs[measurement_model_parameter_index]^2)
  max_measurement_model = dmvnorm(t(matrix(0,meas_dim)), t(matrix(0,meas_dim)), R, log=TRUE)
  
  current_log_llhd_estimate = current_log_llhd_estimate + T*max_measurement_model
  
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
  
  sample[,1] = inputs[index_of_input_in_param_order]
  
  current_inputs_for_simulator = inputs
  
  for (i in 2:num_mcmc)
  {
    if (i%%1==0)
    {
      print(sprintf("Current iteration: %i", i))
      #write.matrix(sample[,1:(i-1)],file=paste("paper_sr02longest_sample_enkf_",num_internal_particles,"_",obs_noise_var,"_",experiment_number,".txt",sep=""))
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
      
      if (length(measurement_model_parameter_index)==1)
        R = diag(proposed_inputs_for_simulator[measurement_model_parameter_index]^2,meas_dim)
      else
        R = diag(proposed_inputs_for_simulator[measurement_model_parameter_index]^2)
      proposed_max_measurement_model = dmvnorm(t(matrix(0,meas_dim)), t(matrix(0,meas_dim)), R, log=TRUE)

      termination_criterion = logu[b,i] - log_accept_prior + T*(-proposed_max_measurement_model) + current_log_llhd_estimate
      
      result = enkf(y,num_particles,proposed_inputs_for_simulator,initial_step_simulator,step_simulator,H,measurement_model_parameter_index,dim_x,resample_threshold,termination_criterion)
      proposed_log_llhd_estimate = result[[1]]
      proposed_data = result[[2]]
      proposed_states = result[[3]]
      
      if (!is.infinite(proposed_log_llhd_estimate))
      {
        if (proposed_log_llhd_estimate>termination_criterion)
        {
          sample[,i] = proposed_parameter
          current_inputs_for_simulator = proposed_inputs_for_simulator
          
          print("Accepted")
          
          if (length(measurement_model_parameter_index)==1)
            R = diag(current_inputs_for_simulator[measurement_model_parameter_index]^2,meas_dim)
          else
            R = diag(current_inputs_for_simulator[measurement_model_parameter_index]^2)
          current_max_measurement_model = dmvnorm(t(matrix(0,meas_dim)), t(matrix(0,meas_dim)), R, log=TRUE)

          current_log_llhd_estimate = proposed_log_llhd_estimate + T*current_max_measurement_model
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
