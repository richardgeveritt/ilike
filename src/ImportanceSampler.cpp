#include <Rcpp.h>
using namespace Rcpp;

#include "UsefulFunctions.h"


// [[Rcpp::export]]
List DoImportanceSamplerCpp(unsigned int NumberOfPoints,
                                     const List &Model,
                                     List Algorithm,
                                     unsigned int MaxVectorSize)
{

  try
  {
    // Work out the number of batches to simulate.
    // We will store the points and the auxiliary variables in files to avoid memory problems.
    NumericVector Inputs = Model["Inputs"];
    IntegerVector ParameterIndex = Model["ParameterIndex"];
    NumericVector Data = Model["Data"];
    //unsigned int InputsDimension = Inputs.length();
    //unsigned int ParameterDimension = Model["ParameterDimension"];
    unsigned int AuxiliaryVariablesDimension = Algorithm["AuxiliaryVariablesDimension"];
    unsigned int TotalDimension = Inputs.length() + AuxiliaryVariablesDimension;
    unsigned int NumberOfBatchesMinusOne = floor((NumberOfPoints*TotalDimension-1)/MaxVectorSize);

    // Simulate batch of parameters and associated aux variables.
    // Save all in file if multiple batches.
    if (NumberOfBatchesMinusOne==0)
    {
      // Everything can be done in one batch, so things are more straightforward.

      // Do the simulation.
      SEXP SimulateProposalSEXP = Algorithm["SimulateProposal"];
      NumericMatrix ProposedPoints = MakeSimulateProposalPtrFromSEXP(SimulateProposalSEXP)(NumberOfPoints);
      NumericMatrix ProposedInputs(NumberOfPoints,Inputs.size());
      for (unsigned int i=0; i<NumberOfPoints; ++i)
      {
        NumericVector ProposedInputsRow = Inputs;
        ProposedInputsRow[ParameterIndex] = NumericVector(ProposedPoints(i,_));
        ProposedInputs(i,_) = ProposedInputsRow;
      }

      List LikelihoodEstimator = Algorithm["LikelihoodEstimator"];

      bool EvaluateLogLikelihoodIsSet = LikelihoodEstimator["EvaluateLogLikelihoodIsSet"];

      List ProposedAuxiliaryVariables;

      if (EvaluateLogLikelihoodIsSet==FALSE)
      {
        SEXP SimulateAuxiliaryVariablesSEXP = LikelihoodEstimator["SimulateAuxiliaryVariables"];
        SimulateAuxiliaryVariablesPtr SimulateAuxiliaryVariables = MakeSimulateAuxiliaryVariablesPtrFromSEXP(SimulateAuxiliaryVariablesSEXP);

        for (unsigned int i=0; i<NumberOfPoints; ++i)
        {
          ProposedAuxiliaryVariables.push_back(SimulateAuxiliaryVariables(ProposedInputs(i,_),Data));
        }

        // Now configure the likelihood estimator using all of the simulations, if needed.
        SEXP SetUpLikelihoodEstimatorSEXP = LikelihoodEstimator["SetUpLikelihoodEstimator"];
        SetUpLikelihoodEstimatorPtr SetUpLikelihoodEstimator = MakeSetUpLikelihoodEstimatorPtrFromSEXP(SetUpLikelihoodEstimatorSEXP);

        XPtr<EvaluateLogLikelihoodPtr> EvaluateLogLikelihood = SetUpLikelihoodEstimator(ProposedPoints,ProposedAuxiliaryVariables);
        SEXP EvaluateLogLikelihoodSEXP = SEXP(EvaluateLogLikelihood);
        LikelihoodEstimator["EvaluateLogLikelihood"] = EvaluateLogLikelihoodSEXP;
        Algorithm["LikelihoodEstimator"] = LikelihoodEstimator;
      }

      // Calculate weights.
      SEXP EvaluateLogLikelihoodSEXP = LikelihoodEstimator["EvaluateLogLikelihood"];
      EvaluateLogLikelihoodPtr EvaluateLogLikelihood = MakeEvaluateLogLikelihoodPtrFromSEXP(EvaluateLogLikelihoodSEXP);

      NumericVector LogWeights(NumberOfPoints);
      bool PriorIsProposal = Algorithm["PriorIsProposal"];
      if (PriorIsProposal==TRUE)
      {
        for (unsigned int i=0; i<NumberOfPoints; ++i)
        {
          LogWeights[i] = EvaluateLogLikelihood(ProposedInputs(i,_),Data);
        }
      }
      else
      {
        SEXP EvaluateLogPriorSEXP = Model["EvaluateLogPrior"];
        EvaluateDistributionPtr EvaluateLogPrior = MakeEvaluateDistributionPtrFromSEXP(EvaluateLogPriorSEXP);

        SEXP EvaluateLogProposalSEXP = Algorithm["EvaluateLogProposal"];
        EvaluateDistributionPtr EvaluateLogProposal = MakeEvaluateDistributionPtrFromSEXP(EvaluateLogProposalSEXP);

        for (unsigned int i=0; i<NumberOfPoints; ++i)
        {
          LogWeights[i] = EvaluateLogLikelihood(NumericVector(ProposedInputs(i,_)),Data) + EvaluateLogPrior(NumericVector(ProposedPoints(i,_))) - EvaluateLogProposal(NumericVector(ProposedPoints(i,_)));
        }
      }

      return List::create(Named("ProposedPoints") = ProposedPoints,
                          Named("ProposedAuxiliaryVariables") = ProposedAuxiliaryVariables,
                          Named("LogWeights") = LogWeights,
                          Named("LogNormalisingConstant") = LogSumExpCpp(LogWeights));
    }
    else
    {
      // need to sort out error handling
      throw 1;
    }
  }
  catch (int e)
  {
    if (e==1)
      Rcout << "Cannot yet deal with multiple batches." << e << '\n';
    return List::create();
  }

}
