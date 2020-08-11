#include <Rcpp.h>
using namespace Rcpp;

#ifndef USEFULFUNCTIONS_H
#define USEFULFUNCTIONS_H

typedef NumericMatrix (*SimulateProposalPtr)(const unsigned int &NumberOfPoints);

SimulateProposalPtr MakeSimulateProposalPtrFromSEXP(const SEXP &SimulateProposalSEXP)
{
  XPtr<SimulateProposalPtr> LlhdSimulateXPtr(SimulateProposalSEXP);
  return *LlhdSimulateXPtr;
}


typedef double (*EvaluateLogLikelihoodPtr)(const NumericVector &Inputs, const NumericVector &y);

EvaluateLogLikelihoodPtr MakeEvaluateLogLikelihoodPtrFromSEXP(const SEXP &EvaluateLogLikelihoodSEXP)
{
  XPtr<EvaluateLogLikelihoodPtr> EvaluateLogLikelihoodXPtr(EvaluateLogLikelihoodSEXP);
  return *EvaluateLogLikelihoodXPtr;
}


typedef List (*SimulateAuxiliaryVariablesPtr)(const NumericVector &Inputs, const NumericVector &y);

SimulateAuxiliaryVariablesPtr MakeSimulateAuxiliaryVariablesPtrFromSEXP(const SEXP &SimulateAuxiliaryVariablesSEXP)
{
  XPtr<SimulateAuxiliaryVariablesPtr> SimulateAuxiliaryVariablesXPtr(SimulateAuxiliaryVariablesSEXP);
  return *SimulateAuxiliaryVariablesXPtr;
}


typedef XPtr<EvaluateLogLikelihoodPtr> (*SetUpLikelihoodEstimatorPtr)(const NumericMatrix &Inputs, const List &AuxiliaryVariables);

SetUpLikelihoodEstimatorPtr MakeSetUpLikelihoodEstimatorPtrFromSEXP(const SEXP &SetUpLikelihoodEstimatorSEXP)
{
  XPtr<SetUpLikelihoodEstimatorPtr> SetUpLikelihoodEstimatorXPtr(SetUpLikelihoodEstimatorSEXP);
  return *SetUpLikelihoodEstimatorXPtr;
}


typedef double (*EvaluateDistributionPtr)(const NumericVector &Inputs);

EvaluateDistributionPtr MakeEvaluateDistributionPtrFromSEXP(const SEXP &EvaluateDistributionSEXP)
{
  XPtr<EvaluateDistributionPtr> EvaluateDistributionXPtr(EvaluateDistributionSEXP);
  return *EvaluateDistributionXPtr;
}

// // From SymbolixAU at https://stackoverflow.com/questions/59055902/find-index-of-all-max-min-values-in-vector-in-rcpp
// IntegerVector which_max(const NumericVector &v)
// {
//   double current_max = v[0];
//   int n = v.size();
//   std::vector< int > res;
//   res.push_back( 0 );
//   int i;
//   
//   for( i = 1; i < n; ++i) {
//     double x = v[i];
//     if( x > current_max ) {
//       res.clear();
//       current_max = x;
//       res.push_back( i );
//     } else if ( x == current_max ) {
//       res.push_back( i );
//     }
//   }
//   Rcpp::IntegerVector iv( res.begin(), res.end() );
//   return iv;
// }

double LogSumExpCpp(const NumericVector &LogWeights)
{
  if (sum(is_na(LogWeights))>0)
  {
    return R_NegInf;
  }
  unsigned int xmax = which_max(LogWeights);
  IntegerVector xmaxv = IntegerVector::create(xmax);
  unsigned int num_weights = LogWeights.length();
  IntegerVector all_but_one = setdiff(IntegerVector(seq(0,num_weights-1)),xmaxv);
  
  NumericVector AllButOneLogWeights = LogWeights[all_but_one];
    
  double result = log1p(sum(exp(AllButOneLogWeights-LogWeights[xmax])))+LogWeights[xmax];
  if (isnan(result))
  {
    return R_NegInf;
  }
  else
  {
    return result;
  }
}

#endif