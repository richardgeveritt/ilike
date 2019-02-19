library(LaplacesDemon)
library(gsubfn)
library("future.apply")
library("future.batchtools")

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

score_u <- function(epsilon,current_error,alpha,num_particles,resample_u,current_u)
{
  which_nonzero = which(current_error<epsilon)
  log_weights = matrix(-Inf,num_particles)
  log_weights[which_nonzero] = -log(num_particles)
  if (all(log_weights==-Inf))
  {
    1-alpha*num_particles
  } else {
    indices = stratified_resample_u(log_weights-logsumexp(log_weights),resample_u)
    current_u = current_u[indices,]
    #nrow(unique(current_u))-alpha*num_particles
    all_unique = sapply(lapply(lapply(as.list(data.frame(current_u)),unique),length),c)
    #print(median(all_unique))
    median(all_unique)-alpha*num_particles
  }
}

bisection_u <- function(range,current_error,alpha,num_particles,resample_u,current_u)
{
  current_bisect_epsilon = range[2]
  bisect_size = (range[2] - range[1])/2
  direction = -1
  
  old_score = score_u(current_bisect_epsilon,current_error,alpha,num_particles,resample_u,current_u)
  
  for (i in 1:100)
  {
    new_bisect_epsilon = current_bisect_epsilon + direction*bisect_size
    bisect_size = bisect_size/2
    the_score = score_u(new_bisect_epsilon,current_error,alpha,num_particles,resample_u,current_u)
    direction = -sign(the_score)
    current_bisect_epsilon = new_bisect_epsilon
    
    if (the_score==0)
    {
      break
    }
    
  }
  
  current_bisect_epsilon
}

update_distance = function(current_distance,distance_type,summary_scaling)
{
  if (distance_type=="none")
  {
    return(current_distance)
  }
  else
  {
    if (distance_type=="euclidean") {
      distance = function(observed,simulated) {
        return(sqrt(sum(sum(((simulated-observed)*summary_scaling)^2))))
      }
    } else {
      # throw error
      browser()
    }
    return(distance)
  }
}

mcmc_prior = function(sample,b,num_within_block,index_of_input_in_param_order,index_of_param_in_input_order,priors_index_in_proposal,proposals_index,priors_index,proposals_simulate,proposals_logevaluate,priors_logevaluate,prop_cov)
{
  proposed_parameter = sample
  log_accept_prior = 0
  proposed_parameter[index_of_param_in_input_order[proposals_index[[b]]]] = proposals_simulate[[b]](conditioned_on_parameter=sample[index_of_param_in_input_order[proposals_index[[b]]]],other_parameters=sample,proposal_covariance=prop_cov)
  for (w in 1:num_within_block[b])
  {
    log_accept_prior = log_accept_prior + priors_logevaluate[[ priors_index_in_proposal[[b]][w] ]](proposed_parameter[ index_of_param_in_input_order[priors_index[[priors_index_in_proposal[[b]][w]]]] ]) - priors_logevaluate[[ priors_index_in_proposal[[b]][w] ]](sample[ index_of_param_in_input_order[priors_index[[priors_index_in_proposal[[b]][w]]]]])
  }
  log_accept_prior = log_accept_prior + proposals_logevaluate[[b]](current_parameter=sample[index_of_param_in_input_order[proposals_index[[b]]]],conditioned_on_parameter=proposed_parameter[ index_of_param_in_input_order[proposals_index[[b]]] ],other_parameters=sample,proposal_covariance=prop_cov) - proposals_logevaluate[[b]](current_parameter=proposed_parameter[ index_of_param_in_input_order[proposals_index[[b]]] ],conditioned_on_parameter=sample[ index_of_param_in_input_order[proposals_index[[b]]] ],other_parameters=sample,proposal_covariance=prop_cov)
  
  return(list(log_accept_prior,proposed_parameter))
}

abc_smc = function(num_particles,num_unique,y,simulator,inputs=matrix(0,0,0),priors_logevaluate,priors_simulate,priors_index=list(),proposals_simulate,proposals_logevaluate,proposals_index=list(),tolerance=Inf,num_abc_sims=1,abc_summary="identity",summary_scaling=1,abc_kernel="uniform",distance="euclidean",store_states=FALSE,is_net_logo=FALSE)
{
  
  alpha = num_unique/num_particles
  # Check dimensions of priors_evaluate,proposals_simulate,proposals_logevaluate,parameter_index
  
  # Case 1 - we have a single prior and single proposal, no inputs and no parameter index
  # In this case we need to set up parameter_index, enums_for_block, inputs, prior/proposal in block, and param count
  # case - we have single prior, and inputs and no param index
  if ( (length(priors_index)==0) || (length(proposals_index)==0) )
  {
    # No parameter_index supplied, so the assumption is that all of the dimensions in the priors, proposals, inputs and simulator refer to a vector in the same order.
    
    if ( (!is.function(priors_logevaluate)) || (!is.function(proposals_simulate)) || (!is.function(proposals_logevaluate)) || (!is.function(priors_simulate)) )
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
    } else if ( ( (!is.matrix(inputs)) && (length(inputs)!=0) ) || ( (is.matrix(inputs)) && (nrow(inputs)!=0) ) ) {
      priors_index = list()
      proposals_index = list()
      priors_index[[1]] = 1:length(inputs)
      proposals_index[[1]] = 1:length(inputs)
    } else {
      inputs = priors_simulate(1)
      priors_index = list()
      priors_index[[1]] = 1:length(inputs)
      proposals_index = list()
      proposals_index[[1]] = 1:length(inputs)
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
  
  if (length(priors_simulate)!=length(priors_logevaluate))
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
  
  if ( ( (!is.matrix(inputs)) && (length(inputs)==0) ) || ( (is.matrix(inputs)) && (nrow(inputs)==0) ) ) {
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
  sim_output = matrix(0,0,0)
  sim_states = matrix(0,0,0)
  res = simulator(inputs)
  sim_output = res[[1]]
  sim_states = res[[2]]
  if (length(sim_output)!=length(y))
  {
    # throw error
    browser()
  }
  
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
    distance_type = distance
  } else {
    distance_type = "none"
  }
  
  sum_y = abc_summary(y)
  
  particles = matrix(0,num_particles,num_params)
  
  if (store_states==TRUE)
    states = matrix(0,num_particles,num_states)  
  else
    states = matrix(0,0,0)
  
  proposed_inputs_for_simulator = t(matrix(rep(inputs,num_particles),length(inputs),num_particles))
  for (j in 1:length(priors_index))
  {
    particles[,index_of_param_in_input_order[priors_index[[j]]]] = priors_simulate[[j]](num_particles)
    proposed_inputs_for_simulator[,priors_index[[j]]] = particles[,index_of_param_in_input_order[priors_index[[j]]]]
  }
  
  current_inputs_for_simulator = proposed_inputs_for_simulator
  
  rows = lapply(1:num_particles,function(i){proposed_inputs_for_simulator[i,]})
  example_output = simulator(proposed_inputs_for_simulator[1,])
  if (is.list(example_output))
  {
    state_output = TRUE
  }
  else
  {
    state_output = FALSE
  }
  
  if (is_net_logo==FALSE)
  {
    sim_data = future_lapply(rows,FUN=function(th){simulator(th)},future.seed = TRUE)
  }
  else
  {
    browser()
  }
  
  if (state_output==TRUE)
  {
    if (store_states==TRUE)
    {
      states = lapply(sim_data, `[[`, 2)
    }
    
    sim_data = lapply(sim_data, `[[`, 1)
  }
  
  sumstats = future_lapply(sim_data,FUN=abc_summary)
  sumstats_vec = sapply(sumstats,c)
  #summary_scaling = 1/apply(sumstats_vec,1,sd,na.rm=TRUE)
  
  distance = update_distance(distance,distance_type,summary_scaling)
  distance_from_y = function(stat)
  {
    distance(stat,sum_y)
  }
  
  dist = future_lapply(sumstats,FUN=distance_from_y)
  dist = sapply(dist,c)
  
  current_epsilon = Inf
  current_w = matrix(-log(num_particles),num_particles)
  
  eps_store = c()
  current_epsilon = max(dist)
  
  while (current_epsilon>tolerance)
  { 
    dist[which(is.na(dist))] = Inf
    old_epsilon = current_epsilon
    
    resample_u = runif(num_particles,0,1)
    
    if (length(which(!is.infinite(dist)))>alpha*num_particles)
      current_epsilon = bisection_u(c(0,current_epsilon),dist,alpha,num_particles,resample_u,particles)
    else
      current_epsilon = max(dist[which(!is.infinite(dist))])
    
    print(sprintf("Current tolerance: %f", current_epsilon))

    if ( (current_epsilon<=tolerance) )
    {
      current_w[which(dist>tolerance)] = -Inf
      eps_store = c(eps_store,tolerance)
      
      indices = stratified_resample_u(current_w-logsumexp(current_w),resample_u)
      sim_data = sim_data[indices]#current_x[indices,,]
      particles = particles[indices,]
      current_inputs_for_simulator = current_inputs_for_simulator[indices,]
      sumstats_vec = sumstats_vec[,indices]
      dist = dist[indices]
      current_w = matrix(-log(num_particles),num_particles)
      if (store_states==TRUE)
        states = states[indices]
      
      write.matrix(eps_store,file="epsilon.txt")
      write.matrix(particles,file=paste("particles",current_epsilon,counter,".txt",sep="_"))
      
      break
    }
    else
    {
      current_w[which(dist>current_epsilon)] = -Inf
      eps_store = c(eps_store,current_epsilon)
      
      write.matrix(eps_store,file="epsilon.txt")
    }
    
    indices = stratified_resample_u(current_w-logsumexp(current_w),resample_u)
    sim_data = sim_data[indices]#current_x[indices,,]
    particles = particles[indices,]
    current_inputs_for_simulator = current_inputs_for_simulator[indices,]
    sumstats_vec = sumstats_vec[,indices]
    dist = dist[indices]
    
    if (store_states==TRUE)
      states = states[indices]
    current_w = matrix(-log(num_particles),num_particles)
    
    logu = matrix(log(runif(num_particles*num_blocks)),num_blocks,num_particles)
    
    update_order = sample(1:num_blocks,num_blocks)
    for (a in 1:num_blocks)
    {
      b = update_order[a]
      
      if (length(proposals_index[[b]])==num_params)
      {
        if (is.matrix(particles))
          prop_cov = cov(particles)
        else
          prop_cov = var(particles)
      }
      else
      {
        current_indices = index_of_param_in_input_order[proposals_index[[b]]]
        conditioned_on_indices = setdiff(1:num_params,current_indices)
        Z = as.data.frame(particles[,conditioned_on_indices])
        fit = lm(particles[,current_indices]~.,data=Z)
        if (is.matrix(particles[,current_indices]))
          prop_cov = cov(residuals(fit))
        else
          prop_cov = var(residuals(fit))
      }
      
      particle_list = as.list(data.frame(t(particles)))
      prior_output = future_lapply(particle_list,FUN=mcmc_prior,b,num_within_block,index_of_input_in_param_order,index_of_param_in_input_order,priors_index_in_proposal,proposals_index,priors_index,proposals_simulate,proposals_logevaluate,priors_logevaluate,prop_cov,future.seed = TRUE)

      proposed_inputs_for_simulator = current_inputs_for_simulator
      proposed_inputs_for_simulator[,index_of_input_in_param_order] = particles
      proposed_block = t(sapply(lapply(prior_output, `[[`, 2),cbind)[index_of_param_in_input_order[proposals_index[[b]]],])
      proposed_inputs_for_simulator[,proposals_index[[b]]] = proposed_block
      proposed_parameter = particles
      proposed_parameter[,index_of_param_in_input_order[proposals_index[[b]]]] = proposed_block
      
      log_accept_prior = sapply(lapply(prior_output, `[[`, 1),cbind)
      
      which_past_prior = which(logu[b,]<log_accept_prior)
      
      if (length(which_past_prior)>0)
      {
        proposed_inputs_for_simulator = proposed_inputs_for_simulator[which_past_prior,]
        
        rows = lapply(1:nrow(proposed_inputs_for_simulator),function(i){proposed_inputs_for_simulator[i,]})
        
        if (is_net_logo==FALSE)
        {
          proposed_sim_data = future_lapply(rows,function(th){simulator(th)},future.seed = TRUE)
        }
        else
        {
          browser()
        }
        
        if (state_output==TRUE)
        {
          if (store_states==TRUE)
          {
            proposed_sim_states = lapply(proposed_sim_data, `[[`, 2)
          }
          
          proposed_sim_data = lapply(proposed_sim_data, `[[`, 1)
        }
        
        proposed_sumstats = future_lapply(proposed_sim_data,FUN=abc_summary)
        proposed_sumstats_vec = sapply(proposed_sumstats,c)
        proposed_dist = future_lapply(proposed_sumstats,FUN=distance_from_y)
        proposed_dist = sapply(proposed_dist,c)
        which_kept = which(proposed_dist<current_epsilon)
        
        print(sprintf("Accepted: %i", length(which_kept)))
        
        if (length(which_kept)>0)
        {
          
          particles[which_past_prior[which_kept],] = proposed_parameter[which_past_prior[which_kept],]
          current_inputs_for_simulator[which_past_prior[which_kept],] = proposed_inputs_for_simulator[which_kept,]
          sim_data[which_past_prior[which_kept]] = proposed_sim_data[which_kept]
          sumstats_vec[,which_past_prior[which_kept]] = proposed_sumstats_vec[,which_kept]
          dist[which_past_prior[which_kept]] = proposed_dist[which_kept]
          if (store_states==TRUE)
            states[which_past_prior[which_kept]] = proposed_sim_states[which_kept]
          
        }
      }
    }
    
    write.matrix(particles,file=paste("particles",current_epsilon,counter,".txt",sep="_"))
    
  }
  
  if (store_states==TRUE)
    states = t(sapply(states,c))
  sim_data = t(sapply(sim_data,c))
  
  if (is_net_logo==TRUE)
  {
    stopCluster(cl)
  }
  
  return(list(particles,sim_data,states,sumstats_vec,dist))
}