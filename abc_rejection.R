library("future.apply")
library("future.batchtools")
library(abctools)

abc_rejection = function(num_points,proportion,y,simulator,inputs=matrix(0,0,0),priors_simulate,priors_index=list(),num_abc_sims=1,abc_summary="identity",abc_kernel="uniform",distance="euclidean",store_states=FALSE,is_net_logo=FALSE)
{
  if (length(priors_index)==0)
  {
    # No parameter_index supplied, so the assumption is that all of the dimensions in the priors, proposals, inputs and simulator refer to a vector in the same order.
    
    if (!is.function(priors_simulate))
    {
      # throw error
    }
    
    if (length(priors_index)>1) {
      # throw error
    }
    
    if (length(priors_index)==1) {
      if (length(inputs)!=length(priors_index[[1]])) {
        # throw error
      }
      else {
        proposals_index = priors_index
      }
    } else if (length(inputs)!=0) {
      priors_index = list()
      priors_index[[1]] = 1:length(inputs)
    } else {
      inputs = priors_simulate(1)
      priors_index = list()
      priors_index[[1]] = 1:length(inputs)
    }
  }
  
  if (!is.list(priors_index))
  {
    prior_input2 = priors_index
    priors_index = list()
    priors_index[[1]] = prior_input2
  }
  
  if (is.function(priors_simulate))
  {
    prior_input = priors_simulate
    priors_simulate = list()
    priors_simulate[[1]] = prior_input
    
    if ( is.list(priors_index) && (length(priors_index)!=1) )
    {
      # throw error
    }
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
    }
  }
  
  proposed_inputs_for_simulator = t(matrix(rep(inputs,num_points),length(inputs),num_points))
  for (j in 1:length(priors_index))
  {
    prior_sims = priors_simulate[[j]](num_points)
    proposed_inputs_for_simulator[,priors_index[[j]]] = prior_sims
  }
  
  rows = lapply(1:num_points,function(i){proposed_inputs_for_simulator[i,]})
  example_output = simulator(proposed_inputs_for_simulator[1,])
  if (is.list(example_output))
  {
    state_output = TRUE
  }
  else
  {
    state_output = FALSE
  }
  
  sim_data = future_lapply(rows,FUN=function(th){simulator(th)},future.seed = TRUE)
  
  if (is.character(abc_summary))
  {
    if (abc_summary=="identity") {
      abc_summary = function(data){return(data)}
    }
    else {
      # throw error
    }
  } else {
    if (!is.function(abc_summary))
    {
      # throw error
    }
  }
  
  sum_y = abc_summary(y)

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
  summary_scaling = 1/apply(sumstats_vec,1,sd,na.rm=TRUE)
  
  if (is.character(abc_kernel))
  {
    if (!(abc_kernel=="uniform")) {
      # throw error
    }
  } else {
    if (!is.function(abc_kernel))
    {
      # throw error
    }
  }
  
  if (is.character(distance))
  {
    if (distance=="euclidean") {
      distance = function(observed,simulated) {
        return(sqrt(sum(sum(((simulated-observed)*summary_scaling)^2))))
      }
    } else {
      # throw error
    }
  } else {
    if (!is.function(distance))
    {
      # throw error
    }
  }
  
  distance_from_y = function(stat)
  {
    distance(stat,sum_y)
  }
  
  dist = future_lapply(sumstats,FUN=distance_from_y)
  dist = sapply(dist,c)
  cutoff = quantile(dist,proportion,na.rm=TRUE)
  which_kept = (dist<cutoff)
  post = prior_sims[which_kept,]
  
  if (store_states==TRUE)
    states = t(sapply(states,c))
  else
    states = matrix(0,0,0)

  sim_data = t(sapply(sim_data,c))
  
  # Perform F&P dim reduction
  fandp_info = saABC(prior_sims,sim_data,plot=FALSE)
  fandp_sumstats_vec = sim_data%*%t(fandp_info$B)
  fandp_summary_scaling = 1/apply(fandp_sumstats_vec,2,sd,na.rm=TRUE)
  
  return(list(post,cutoff,which_kept,prior_sims,sim_data,states,sumstats_vec,dist,summary_scaling,index_of_input_in_param_order,fandp_info,fandp_sumstats_vec,fandp_summary_scaling))
}