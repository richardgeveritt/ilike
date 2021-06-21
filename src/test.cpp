#include <RcppArmadillo.h>
using namespace Rcpp;

#include <RcppParallel.h>

#include "parameters.h"
#include "function_pointers.h"



// class Tester {
//
// public:
//
//   Tester(double input) {stored = input;}
//
//   double stored;
//
//   double get() { return(stored/2.0); }
// };


inline double my_sqrt(const Parameters &badger)
{
  arma::colvec thing = badger["theta"];
  sleep(5);
  //Tester thing(20.0);

  return badger["theta"][2]/2.0;
}

class thing {

public:

  thing() {}
  ~thing() {}

  double operator()(const Parameters &inputs) { double output = sum(inputs["theta"]);
    return output; }

};


class MyWorker : public RcppParallel::Worker {

public:

  MyWorker(const std::vector<Parameters> &input_in,
           const std::vector<double> &dummy_input_in,
           const EvaluateLogDistributionPtr &evaluate_log_prior_in,
           const SimulateDistributionPtr &simulate_prior_in);

  void operator()(std::size_t begin, std::size_t end) {
    std::transform(this->input.begin() + begin,
                   this->input.begin() + end,
                   this->output.begin() + begin,
                   this->evaluate_log_prior);
    // std::transform(this->dummy_input.begin() + begin,
    //                this->dummy_input.begin() + end,
    //                this->simulate_output.begin() + begin,
    //                this->simulate_prior);
  }

  std::vector<double> output;

  std::vector<Parameters> simulate_output;

  thing my_thing;

private:

  std::vector<Parameters> input;

  std::vector<double> dummy_input;

  EvaluateLogDistributionPtr evaluate_log_prior;

  SimulateDistributionPtr simulate_prior;

};

MyWorker::MyWorker(const std::vector<Parameters> &input_in,
                   const std::vector<double> &dummy_input_in,
                   const EvaluateLogDistributionPtr &evaluate_log_prior_in,
                   const SimulateDistributionPtr &simulate_prior_in)
{
  this->input = input_in;
  this->dummy_input = dummy_input_in;
  this->output = std::vector<double>(this->input.size());
  this->simulate_output = std::vector<Parameters>(this->dummy_input.size());
  this->evaluate_log_prior = evaluate_log_prior_in;
  this->simulate_prior = simulate_prior_in;
}

// [[Rcpp::export]]
double a_test(const List &model)
{
  SEXP evaluate_log_prior_SEXP = model["evaluate_log_prior"];
  EvaluateLogDistributionPtr evaluate_log_prior = load_evaluate_log_distribution(evaluate_log_prior_SEXP);

  SEXP simulate_prior_SEXP = model["simulate_prior"];
  SimulateDistributionPtr simulate_prior = load_simulate_distribution(simulate_prior_SEXP);

  Parameters test_params = simulate_prior(1.0);
  Rcout << test_params << std::endl;
  // arma::colvec my_vec(10, arma::fill::zeros);
  // //Rcout << "The value is:" << std::endl << my_vec << std::endl;
  // my_vec[1] = 2.0;
  // test_params["theta"] = my_vec;
  // test_params["theta"][2] = 5.0;
  // //Rcout << "The value is:" << std::endl << test_params["theta"] << std::endl;

  //std::vector<Parameters> input_vec(2);
  //input_vec[0] = test_params;
  //input_vec[1] = test_params;

  //double testing = 4.0;//my_sqrt(test_params);

  std::vector<Parameters> input_vec(10);
  input_vec[0] = test_params;
  input_vec[1] = test_params;
  input_vec[2] = test_params;
  input_vec[3] = test_params;
  input_vec[4] = test_params;
  input_vec[5] = test_params;
  input_vec[6] = test_params;
  input_vec[7] = test_params;
  input_vec[8] = test_params;
  input_vec[9] = test_params;

  std::vector<double> dummy_input_vec(10);
  dummy_input_vec[0] = 0.0;
  dummy_input_vec[1] = 0.0;
  dummy_input_vec[2] = 0.0;
  dummy_input_vec[3] = 0.0;
  dummy_input_vec[4] = 0.0;
  dummy_input_vec[5] = 0.0;
  dummy_input_vec[6] = 0.0;
  dummy_input_vec[7] = 0.0;
  dummy_input_vec[8] = 0.0;
  dummy_input_vec[9] = 0.0;

  MyWorker a_worker(input_vec, dummy_input_vec, evaluate_log_prior, simulate_prior);

  parallelFor(0, dummy_input_vec.size(), a_worker);
  //double answer = 0.0;
  //for (unsigned int i=0; i<10; ++i)
  //  answer = evaluate_log_prior(test_params);
  return a_worker.output[0];
}
