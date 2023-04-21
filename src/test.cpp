#include <RcppArmadillo.h>
using namespace Rcpp;

#include <RcppParallel.h>

/*
#include "parameters.h"
#include "function_pointers.h"
#include "distributions.h"
#include "importance_sampler.h"
#include "likelihood_estimator_output.h"
#include "smc_output.h"
*/


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


// inline double my_sqrt(const Parameters &badger)
// {
//   arma::colvec thing = badger["theta"];
//   sleep(5);
//   //Tester thing(20.0);
//
//   return badger["theta"][2]/2.0;
// }
//
// class thing {
//
// public:
//
//   thing() {}
//   ~thing() {}
//
//   double operator()(const Parameters &inputs) { double output = sum(inputs["theta"]);
//     return output; }
//
// };
//
//

/*
class MyWorker : public RcppParallel::Worker {

public:

  MyWorker(uint64_t seed_in,
           const std::vector<Parameters> &input_in,
           const std::vector<double> &dummy_input_in,
           const EvaluateLogDistributionPtr &evaluate_log_prior_in,
           const SimulateDistributionPtr &simulate_prior_in);

  void operator()(std::size_t begin, std::size_t end) {
    // std::transform(this->input.begin() + begin,
    //                this->input.begin() + end,
    //                this->output.begin() + begin,
    //                this->evaluate_log_prior);
    RandomNumberGenerator local_rng(rng);
    local_rng.seed(seed,end);
    for (std::size_t i = begin; i < end; ++i)
    {
      this->simulate_output[i] = this->simulate_prior(this->rng);
    }
  }

  std::vector<double> output;

  std::vector<Parameters> simulate_output;

private:

  uint64_t seed;

  std::vector<Parameters> input;

  std::vector<double> dummy_input;

  EvaluateLogDistributionPtr evaluate_log_prior;

  SimulateDistributionPtr simulate_prior;

  RandomNumberGenerator rng;

};

MyWorker::MyWorker(uint64_t seed_in,
                   const std::vector<Parameters> &input_in,
                   const std::vector<double> &dummy_input_in,
                   const EvaluateLogDistributionPtr &evaluate_log_prior_in,
                   const SimulateDistributionPtr &simulate_prior_in)
{
  this->seed = seed_in;
  this->input = input_in;
  //this-rng = dqrng::random_64bit_wrapper<dqrng::xoshiro256plus>();
  this->dummy_input = dummy_input_in;
  this->output = std::vector<double>(this->input.size());
  this->simulate_output = std::vector<Parameters>(this->dummy_input.size());
  this->evaluate_log_prior = evaluate_log_prior_in;
  this->simulate_prior = simulate_prior_in;
}
*/

// [[Rcpp::export]]
double a_test(const List &model)
{
  //SEXP evaluate_log_prior_SEXP = model["evaluate_log_prior"];
  //EvaluateLogDistributionPtr evaluate_log_prior = load_evaluate_log_distribution(evaluate_log_prior_SEXP);

  /*
  SEXP simulate_prior_SEXP = model["simulate_prior"];
  SimulateIndependentProposalPtr simulate_prior_func = load_simulate_independent_proposal(simulate_prior_SEXP);

  SEXP evaluate_log_likelihood_SEXP = model["evaluate_log_likelihood"];
  EvaluateLogLikelihoodPtr evaluate_log_likelihood_func = load_evaluate_log_likelihood(evaluate_log_likelihood_SEXP);

  SEXP data_SEXP = model["data"];
  Data data = load_data(data_SEXP);
  Data* data_pointer = &data;
  Rcout << data << std::endl;

  RandomNumberGenerator rng;
  size_t seed = rdtsc();

  bool parallel = FALSE;
  bool smcfixed_flag = TRUE;
  size_t grain_size = 1;

  size_t number_of_particles = 1000000;

  ImportanceSampler is(&rng,
                       &seed,
                       data_pointer,
                       number_of_particles,
                       evaluate_log_likelihood_func,
                       simulate_prior_func,
                       smcfixed_flag,
                       parallel,
                       grain_size);

  clock_t start, end;
  int max = 0;
  start = clock();


  SMCOutput* output = is.run();

  end = clock();
  printf ("time: %0.8f sec, max = %d\n",
          ((float) end - start)/CLOCKS_PER_SEC, max);


  Rcout << output->log_likelihood << std::endl;
  delete output;
  */

  //Parameters input;
  //dqrng::random_64bit_wrapper<dqrng::xoshiro256plus> rng = dqrng::random_64bit_wrapper<dqrng::xoshiro256plus>();
  //rng.seed(1,3);
  //Parameters test_params = simulate_prior(rng);
  //Rcout << test_params << std::endl;
  // Parameters test_params2 = simulate_prior(rng);
  // Rcout << test_params2 << std::endl;
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

  //std::vector<Parameters> input_vec(10);
  // input_vec[0] = test_params;
  // input_vec[1] = test_params;
  // input_vec[2] = test_params;
  // input_vec[3] = test_params;
  // input_vec[4] = test_params;
  // input_vec[5] = test_params;
  // input_vec[6] = test_params;
  // input_vec[7] = test_params;
  // input_vec[8] = test_params;
  // input_vec[9] = test_params;
  //
  //std::vector<double> dummy_input_vec(10);
  // dummy_input_vec[0] = 0.0;
  // dummy_input_vec[1] = 0.0;
  // dummy_input_vec[2] = 0.0;
  // dummy_input_vec[3] = 0.0;
  // dummy_input_vec[4] = 0.0;
  // dummy_input_vec[5] = 0.0;
  // dummy_input_vec[6] = 0.0;
  // dummy_input_vec[7] = 0.0;
  // dummy_input_vec[8] = 0.0;
  // dummy_input_vec[9] = 0.0;
  //
  //MyWorker a_worker(42, input_vec, dummy_input_vec, evaluate_log_prior, simulate_prior);
  //parallelFor(0, dummy_input_vec.size(), a_worker);
  // //double answer = 0.0;
  // //for (unsigned int i=0; i<10; ++i)
  // //  answer = evaluate_log_prior(test_params);
  // return a_worker.output[0];

  // for (unsigned int i=0; i<10; ++i)
  //   Rcout << a_worker.simulate_output[i] << std::endl;
  //
  // Rcout << std::endl << std::endl;
  //
  // parallelFor(0, dummy_input_vec.size(), a_worker);
  //
  // for (unsigned int i=0; i<10; ++i)
  //   Rcout << a_worker.simulate_output[i] << std::endl;
  return 1.0;
}
