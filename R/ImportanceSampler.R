#require(future.apply)
#require(Rcpp)
Rcpp::sourceCpp('src/ImportanceSampler.cpp')

source("R/UsefulFunctions.R")
source("R/ABCLikelihood.R")



CheckModelDimensionConsistentWithSimulation = function(StoredParameterDimension,SampleParameter)
{
  SampleParameterDimension = length(SampleParameter)
  if (is.null(StoredParameterDimension)) {
    StoredParameterDimension = SampleParameterDimension
  }
  else {
    if (StoredParameterDimension != SampleParameterDimension) {
      stop("Simulated parameter is inconsistent with information in the specified model about the parameter dimension.")
    }
  }

  return(StoredParameterDimension)
}

CheckInputs = function(Model,Algorithm)
{
  # Check that the "inputs" and parameter index are consistent.
  if ( (is.null(Model$Inputs) ) && (!is.null(Model$ParameterIndex) ) ) {
    # If inputs are not set, simulate them from the proposal.
    Model$Inputs = Algorithm$SimulateProposal(1)
  } else if ( (!is.null(Model$Inputs) ) && (is.null(Model$ParameterIndex) ) ) {
    # If parameter index is not set, simply index all of the inputs.
    Model$ParameterIndex = 1:length(Model$Inputs)
  } else if ( (!is.null(Model$Inputs) ) && (!is.null(Model$ParameterIndex) ) ) {
    if (length(Model$ParameterIndex)>length(Model$Inputs) ) {
      stop("Model$ParameterIndex should be <= the length of Model$Inputs.")
    }
  } else {
    Model$Inputs = Algorithm$SimulateProposal(1)
    Model$ParameterIndex = 1:length(Model$Inputs)
  }

  return(Model)
}


CheckIS = function(Model, Algorithm)
{
  # Check that any information in the Model about the dimension of the parameter is consistent.
  if ( (is.null(Model$ParameterDimension) ) && (!is.null(Model$ParameterIndex) ) ) {
    StoredParameterDimension = length(Model$ParameterIndex)
  } else if ( (!is.null(Model$ParameterDimension) ) && (is.null(Model$ParameterIndex) ) ) {
    StoredParameterDimension = Model$ParameterDimension
  } else if ( (!is.null(Model$ParameterDimension) ) && (!is.null(Model$ParameterIndex) ) ) {
    if (length(Model$ParameterIndex)==Model$ParameterDimension) {
      StoredParameterDimension = Model$ParameterDimension
    } else {
      stop("Model$ParameterDimension is not the same as the length of Model$ParameterIndex.")
    }
  } else {
    if (!is.null(Model$Inputs) )
      StoredParameterDimension = length(Model$Inputs)
    else
      StoredParameterDimension = NULL
  }

  # Simulation cases:
  # Neither prior nor proposal simulator defined.
  # Prior simulator defined; proposal not.
  # Proposal simulator defined; prior not.
  # Prior and proposal simulator defined.
  if (is.null(Algorithm$SimulateProposal))
  {
    if (is.null(Model$SimulatePrior))
    {
      # Neither prior nor proposal simulator defined.
      stop("To use an importance sampler, either Model$SimulatePrior or Proposal$SimulateProposal need to be specified.")
    }
    else
    {
      # Prior simulator defined; proposal not.
      # In this case, the IS will use the prior as the proposal.
      Algorithm$PriorIsProposal = TRUE
      Algorithm$SimulateProposal = Model$SimulatePrior
      print("No method specified for simulating from the proposal. Using the prior as the proposal.")

    }
  }
  else
  {
    if (!is.null(Model$SimulatePrior))
    {
      # Both prior and proposal simulators are defined, so we should check that their dimensions match.
      SampleParameter = Model$SimulatePrior(1)
      StoredParameterDimension = CheckModelDimensionConsistentWithSimulation(StoredParameterDimension,SampleParameter)
    }

    # We have a standard importance sampler with a different proposal and prior.
    Algorithm$PriorIsProposal = FALSE
    print("Importance points will be simuulate from the specified proposal.")

  }

  # Check all parameter dimensions match.
  SampleParameter = Algorithm$SimulateProposal(1)
  StoredParameterDimension = CheckModelDimensionConsistentWithSimulation(StoredParameterDimension,SampleParameter)

  Model$ParameterDimension = StoredParameterDimension

  Model = CheckInputs(Model,Algorithm)

  # Check that we can evaluate the prior and the proposal, and that they take as their argument something of dimension length(inputs).
  if (Algorithm$PriorIsProposal==FALSE)
  {
    if (is.null(Algorithm$EvaluateLogProposal))
    {
      stop("If the prior is not used as the proposal, you need to specify Algorithm$EvaluateLogProposal.")
    }
    else
    {
      tryCatch(Algorithm$EvaluateLogProposal(SampleParameter),error = function(e) {stop("Algorithm$EvaluateLogProposal generates an error when used on a vector of the dimension of the parameter.")})
    }

    if (is.null(Model$EvaluateLogPrior))
    {
      stop("If the prior is not used as the proposal, you need to specify Algorithm$EvaluateLogProposal.")
    }
    else
    {
      tryCatch(Model$EvaluateLogPrior(SampleParameter),error = function(e) {stop("Model$EvaluateLogPrior generates an error when used on a vector of the dimension of the parameter.")})
    }
  }

  # Check that the method for evaluating the likelihood makes sense, and that it takes something of the right dimension.

  if (is.null(Model$Data))
  {
    stop("Model$Data not specified.")
  }

  # Check which methods for evaluating the likelihood are available given the specification.
  # Method 1: analytic likelihood.
  # Method 2: standard ABC
  LikelihoodMethods = matrix(0,2)

  if (!is.null(Model$EvaluateLogLikelihood))
  {
    LikelihoodMethods[1] = 1
  }

  if (!is.null(Model$Simulate))
  {
    LikelihoodMethods[2] = 1
  }

  # Automatically find a method by checking LikelihoodMethods.
  if (is.null(Model$LikelihoodMethod))
  {
    if (LikelihoodMethods[1] == 1) {
      print("Model$EvaluateLogLikelihood is specified, so we will use this to evaluate the likelihood.")
      Model$LikelihoodMethod = "analytic"
    } else if (LikelihoodMethods[2] == 1) {
      print("Model$Simulate is specified, so we will use ABC as an approximate likelihood.")
      Model$LikelihoodMethod = "abc"
    }

  }


  if (!is.null(Model$LikelihoodMethod))
  {
    if (Model$LikelihoodMethod=="analytic")
    {
      if (LikelihoodMethods[1] == 1)
      {
        tryCatch(Model$EvaluateLogLikelihood(Model$Inputs,Model$Data),error = function(e) {stop("Model$EvaluateLogLikelihood generates an error when used on a vector of dimension Model$Inputs.")})
        LikelihoodEstimator = list(EvaluateLogLikelihood = Model$EvaluateLogLikelihood)
      }
      else
      {
        stop("Model$LikelihoodMethod is analytic, but Model$EvaluateLogLikelihood is not defined.")
      }
    }

    if (Model$LikelihoodMethod=="abc")
    {
      if (LikelihoodMethods[2] == 1)
      {
        tryCatch(Model$Simulate(Model$Inputs),error = function(e) {stop("Model$Simulate generates an error when used on a vector of dimension Model$Inputs.")})

        if (is.null(Model$LikelihoodOptions))
        {
          # Need to specifiy function set up likelihood estimate. Info needs to be stored in a class-like object.
          LikelihoodEstimator = list(SimulateAuxiliaryVariables = Model$Simulate,
                                   SetUpLikelihoodEstimator = ABC$SetUpLikelihoodEstimator)
        }
        else
        {
          stop("Not written yet. Should be checking for: M; epsilon; distance.")
        }

        # Estimator itself needs to take aux variables as an arg, then do the calculation using a method specified in a different file, all of which should be included at the top of this file. Specified after we have sampled all of the thetas, to set thresholds, etc.

      }
      else
      {
        stop("Model$LikelihoodMethod is abc, but Model$Simulate is not defined.")
      }

      # Also should check ABC info.
    }

  }
  else
  {
    stop("No method for evaluating or estimating the likelihood is specified.")
  }

  Algorithm$LikelihoodEstimator = LikelihoodEstimator

  if (Model$LikelihoodMethod!="analytic")
  {
    # Test the generation of the auxiliary variables and store their dimension.
    AuxiliaryVariables = tryCatch(Algorithm$LikelihoodEstimator$SimulateAuxiliaryVariables(Model$Inputs,Model$Data),error = function(e) {stop("Algorithm$LikelihoodEstimator$SimulateAuxiliaryVariables generates an error when used on a vector of dimension Model$Inputs.")})

    Algorithm$AuxiliaryVariablesDimension = length(unlist(AuxiliaryVariables))
  }
  else
  {
    Algorithm$AuxiliaryVariablesDimension = 0
  }


  # Need to check that Algorithm$SimulateProposal(2) gives something of the right dimension.

  # Need to check that all of the Algorithm$ parts exist.

  list(Model=Model,Algorithm=Algorithm)

}


SimulateBatch = function(BatchNumber,NumBatchPoints,Model,Algorithm)
{
  ProposedPoints = Algorithm$SimulateProposal(NumBatchPoints)
  ProposedInputs = t(matrix(rep(Model$Inputs,NumBatchPoints),length(Model$Inputs),NumBatchPoints))
  ProposedInputs[,Model$ParameterIndex] = ProposedPoints
  ProposedInputs = lapply(1:NumBatchPoints,function(i){ProposedInputs[i,]})
  if (!is.null(Algorithm$LikelihoodEstimator))
  {
    ProposedAuxiliaryVariables = future_lapply(ProposedInputs,FUN=function(Input){Algorithm$LikelihoodEstimator$SimulateAuxiliaryVariables(Model$Inputs,Model$Data)},future.seed = TRUE)
  }
}

#' Importance sampler.
#'
#' @param NumberOfPoints The number of importance points.
#' @param Model A model.
#' @param Algorithm Algorithm details.
#' @return Importance points.
#' @examples
#' ImportanceSampler(10000,MyModel,MyAlgorithm)
#' @export
ImportanceSampler = function(NumberOfPoints,
                             Model,
                             Algorithm,
                             MaxVectorSize = 1e+32)
{
  # Check that inputs make sense.
  # Need also to:
  # - get dimensions of inputs/parameters
  # - get dimensions of aux variables
  # - make sure indexing of inputs/parameters is stored if necessary
  # - make sure prior is in standard format and store it
  Output = CheckIS(Model,Algorithm)
  Model = Output$Model
  Algorithm = Output$Algorithm

  # Work out the number of batches to simulate.
  # We will store the points and the auxiliary variables in files to avoid memory problems.
  TotalDimension = length(Model$Inputs) + Algorithm$AuxiliaryVariablesDimension
  NumberOfBatchesMinusOne = floor((NumberOfPoints*TotalDimension-1)/MaxVectorSize)

  # Simulate batch of parameters and associated aux variables.
  # Save all in file if multiple batches
  # Assume proposal is easy to simulate, and use user-provided function for simulating n times, rather than parallelising.
  # Assume it is expensive to generate the auxiliary variables, and parallelise over this.
  if (NumberOfBatchesMinusOne==0)
  {
    # Everything can be done in one batch, so things are more straightforward.

    # Do the simulation.
    ProposedPoints = Algorithm$SimulateProposal(NumberOfPoints)
    ProposedInputs = t(matrix(rep(Model$Inputs,NumberOfPoints),length(Model$Inputs),NumberOfPoints))
    ProposedInputs[,Model$ParameterIndex] = ProposedPoints
    ProposedInputs = lapply(1:NumberOfPoints,function(i){ProposedInputs[i,]})
    if (is.null(Algorithm$LikelihoodEstimator$EvaluateLogLikelihood))
    {
      ProposedAuxiliaryVariables = future_lapply(ProposedInputs,FUN=function(Input){Algorithm$LikelihoodEstimator$SimulateAuxiliaryVariables(Input,Model$Data)},future.seed = TRUE)

      # Now configure the likelihood estimator using all of the simulations, if needed.
      Algorithm$LikelihoodEstimator$EvaluateLogLikelihood = tryCatch(Algorithm$LikelihoodEstimator$SetUpLikelihoodEstimator(ProposedPoints,ProposedAuxiliaryVariables),error = function(e) {stop("Algorithm$LikelihoodEstimator$SetUpLikelihoodEstimator throws an error when used on the proposed points.")})
    }
    else
    {
      ProposedAuxiliaryVariables = NULL
    }

    # Calculate weights.
    if (Algorithm$PriorIsProposal==TRUE)
    {
      LogWeights = unlist(future_lapply(ProposedInputs,function(i){ Algorithm$LikelihoodEstimator$EvaluateLogLikelihood(i,Model$Data) }))
    }
    else
    {
      LogWeights = unlist(future_lapply(ProposedInputs,function(p){ p = i[Model$ParameterIndex]; Model$EvaluateLogPrior(p) - Algorithm$EvaluateLogProposal(p) + Algorithm$LikelihoodEstimator$EvaluateLogLikelihood(i,Model$Data) }))
    }

    Results = list(ProposedPoints = ProposedPoints,
                   ProposedAuxiliaryVariables = ProposedAuxiliaryVariables,
                   LogWeights = LogWeights,
                   LogNormalisingConstant = LogSumExp(LogWeights))

  }
  else
  {
    stop("ImportanceSampler not yet set upt to use multiple batches.")

    # We need to use multiple batches since we don't want to store big vectors in memory.
    MostBatchSizes = floor(NumberOfPoints/NumberOfBatchesMinusOne)
    LastBatchSize = ((NumberOfPoints*TotalDimension)%%MaxVectorSize)/TotalDimension

    # This will store the results to a file.
    lapply(1:NumberOfBatchesMinusOne,FUN=function(i) {SimulateBatch(i,MostBatchSizes,Model,Algorithm); return(NULL)})
    SimulateBatch(NumberOfBatchesMinusOne+1,LastBatchSize,Model,Algorithm)
  }

}


#' Importance sampler using cpp functions.
#'
#' @param NumberOfPoints The number of importance points.
#' @param Model A model.
#' @param Algorithm Algorithm details.
#' @return Importance points.
#' @examples
#' ImportanceSampler(10000,MyModel,MyAlgorithm)
#' @export
ImportanceSamplerCpp = function(NumberOfPoints,
                                Model,
                                Algorithm,
                                MaxVectorSize = 1e+32)
{
  # Check that inputs make sense.
  # Need also to:
  # - get dimensions of inputs/parameters
  # - get dimensions of aux variables
  # - make sure indexing of inputs/parameters is stored if necessary
  # - make sure prior is in standard format and store it
  Output = CheckIS(Model,Algorithm)
  Model = Output$Model
  Algorithm = Output$Algorithm

  DoImportanceSamplerCpp(NumberOfPoints,
                         Model,
                         Algorithm,
                         MaxVectorSize)
}
