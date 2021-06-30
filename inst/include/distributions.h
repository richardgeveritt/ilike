#ifndef SIMULATION
#define SIMULATION

//#include <random>
//#include <boost/random/binomial_distribution.hpp>
//#include <boost/random/mersenne_twister.hpp>
//#include <boost/random/normal_distribution.hpp>
#include <math.h>
#include <dqrng_distribution.h>
#include <dqrng_generator.h>
#include <boost/random/gamma_distribution.hpp>
#include <boost/math/distributions.hpp>
//#include <boost/math/distributions/gamma.hpp>
//#include <boost/math/distributions/normal.hpp>
#include <boost/math/special_functions/gamma.hpp>

//#include <pcg_random.hpp>
//#include <xoshiro.h>
#include <RcppCommon.h>



#define BOOST_DISABLE_ASSERTS 1

//using RandomNumberGenerator = dqrng::random_64bit_wrapper<dqrng::xoshiro256plus>;
typedef dqrng::random_64bit_wrapper<dqrng::xoshiro256plus> RandomNumberGenerator;
//using Binomial = boost::random::binomial_distribution<int>;


using Gamma = boost::math::gamma_distribution<double>;
//using RNG = dqrng::xoshiro256plus;

inline double normal_simulate(RandomNumberGenerator &rng, double mean, double sd)
{
  dqrng::normal_distribution my_normal;
  return my_normal(rng, dqrng::normal_distribution::param_type(mean, sd));
}

inline double normal_logpdf(double x, double mean, double sd)
{
  if (sd<0)
  {
    return NAN;
  }
  if (sd==0)
  {
    if (x==mean)
      return arma::datum::inf;
    else
      return -arma::datum::inf;
  }
  return -log(sd) - 0.5*log(2.0*M_PI) - 0.5*pow((x-mean)/sd,2.0);
}

inline arma::colvec normal_logpdf(const arma::colvec &x, double mean, double sd)
{
  size_t n = x.size();
  arma::colvec result(n);

  if (sd<0)
  {
    for (size_t i = 0; i<n; ++i)
      result[i] = NAN;
    return result;
  }
  if (sd==0)
  {
    for (size_t i = 0; i<n; ++i)
    {
      if (x[i]==mean)
      {
        result[i] = arma::datum::inf;
      }
      else
      {
        result[i] = -arma::datum::inf;
      }
    }
    return result;
  }

  for (size_t i = 0; i<n; ++i)
    result[i] = -log(sd) - 0.5*log(2.0*M_PI) - 0.5*pow((x[i]-mean)/sd,2.0);
  return result;

}

inline double gamma_simulate(RandomNumberGenerator &rng, double shape, double scale)
{
  boost::random::gamma_distribution<double> my_gamma(shape, scale);
  return my_gamma(rng);
}

inline double gamma_logpdf(double x, double shape, double scale)
{
  if ( (shape<=0) || (scale<=0) )
    return double(NAN);
  if (x<0)
    return -arma::datum::inf;
  return -boost::math::lgamma<double>(shape) - shape*log(scale) + (shape-1.0)*log(x) - x/scale;
}

// namespace Rcpp {
//   SEXP wrap(const RNG& rng);
//   RNG& as(SEXP ptr_RNG);
// }

//RCPP_EXPOSED_CLASS(pcg64);

//#include <Rcpp.h>

// namespace Rcpp {
//
//   inline SEXP wrap(RNG& rng) {
//     return Rcpp::XPtr<RNG>(new RNG(&rng));
//   }
//
//   inline RNG& as(SEXP ptr_RNG) {
//     Rcpp::XPtr<RNG> ptr(ptr_RNG);
//     RNG& rng = *ptr;
//     return rng;
//   }
// }

// #include <array>
// #include <mystdint.h>
// #include <functional>
// #include <algorithm>
// #include <type_traits>
//
// //template<size_t N, int_fast8_t A, int_fast8_t B, int_fast8_t C>
// class MyRNG {
//
// public:
//
//   //using result_type = uint64_t;
//
// private:
//
//   // int_fast8_t A;
//   // int_fast8_t B;
//   // int_fast8_t C;
//   //
//   // std::array<uint64_t, 4> state;
//   //
//   // struct SplitMix {
//   //   SplitMix(const uint64_t& k) : state(k) {}
//   //
//   //   uint64_t operator() () {
//   //     uint64_t z = (state += 0x9e3779b97f4a7c15ULL);
//   //     z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
//   //     z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
//   //     return z ^ (z >> 31);
//   //   }
//   //
//   // private:
//   //   uint64_t state;
//   // };
//   //
//   // uint64_t rotl(const uint64_t x, int k);
//   //
//   // uint64_t next();
//
// public:
//   // inline static constexpr uint64_t min() {return 0.0;};
//   // inline static constexpr uint64_t max() {return UINT64_MAX;};
//
//   MyRNG();
//   virtual ~MyRNG();
//   //xoshiro256plus(uint64_t _seed = 0x85c6ea9eb065ebeeULL);
//
//   // void seed(std::function<uint64_t(void)> rng);
//   //
//   // void seed(uint64_t _seed);
//   //
//   // uint64_t operator() ();
//   //
//   // void jump();
//   // void jump(uint64_t n);
//   // void long_jump();
//   // void long_jump(uint64_t n);
// };
//
// RCPP_EXPOSED_WRAP(MyRNG);
//
// #include <Rcpp.h>
//
//
// // inline uint64_t xoshiro256plus::rotl(const uint64_t x, int k) {
// //   return (x << k) | (x >> (64 - k));
// // }
// //
// // inline uint64_t xoshiro256plus::next() {
// //   size_t N_ = 4;
// //   const uint64_t result = state[0] + state[N_ - 1];
// //
// //   const uint64_t t = state[1] << A;
// //
// //   state[2] ^= state[0];
// //   state[3] ^= state[1];
// //   state[1] ^= state[2];
// //   state[0] ^= state[3];
// //
// //   state[2] ^= t;
// //
// //   state[3] = rotl(state[3], B);
// //
// //   return result;
// // }
//
// inline MyRNG::MyRNG()
// {
//   //seed(_seed);
//   //this->A = 17;
//   //this->B = 45;
//   //this->C = 0;
// }
//
// inline MyRNG::~MyRNG()
// {
//   //seed(_seed);
//   //this->A = 17;
//   //this->B = 45;
//   //this->C = 0;
// }
//
// // inline xoshiro256plus::xoshiro256plus(uint64_t _seed) {
// //   //seed(_seed);
// //   //this->A = 17;
// //   //this->B = 45;
// //   //this->C = 0;
// // }
//
// // inline void xoshiro256plus::seed(std::function<uint64_t(void)> rng) {
// //   std::generate(state.begin(), state.end(), rng);
// // }
// //
// // inline void xoshiro256plus::seed(uint64_t _seed) {
// //   seed(SplitMix(_seed));
// // }
// //
// // inline uint64_t xoshiro256plus::operator() () {
// //   return next();
// // }
// //
// // inline void xoshiro256plus::jump(uint64_t n) {
// //   for( ; n > 0; --n) jump();
// // }
// //
// // inline void xoshiro256plus::long_jump(uint64_t n) {
// //   for( ; n > 0; --n) long_jump();
// // }
// //
// // /* This is xoshiro256+ 1.0, our best and fastest generator for floating-point
// //  numbers. We suggest to use its upper bits for floating-point
// //  generation, as it is slightly faster than xoshiro256**. It passes all
// //  tests we are aware of except for the lowest three bits, which might
// //  fail linearity tests (and just those), so if low linear complexity is
// //  not considered an issue (as it is usually the case) it can be used to
// //  generate 64-bit outputs, too.
// //  We suggest to use a sign test to extract a random Boolean value, and
// //  right shifts to extract subsets of bits.
// //  The state must be seeded so that it is not everywhere zero. If you have
// //  a 64-bit seed, we suggest to seed a splitmix64 generator and use its
// //  output to fill s. */
// //
// // /* This is the jump function for the generator. It is equivalent
// //  to 2^128 calls to next(); it can be used to generate 2^128
// //  non-overlapping subsequences for parallel computations. */
// // inline void xoshiro256plus::jump() {
// //   static const uint64_t JUMP[] = { 0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c };
// //
// //   uint64_t s0 = 0;
// //   uint64_t s1 = 0;
// //   uint64_t s2 = 0;
// //   uint64_t s3 = 0;
// //   for(unsigned int i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
// //     for(unsigned int b = 0; b < 64; b++) {
// //       if (JUMP[i] & uint64_t(1) << b) {
// //         s0 ^= state[0];
// //         s1 ^= state[1];
// //         s2 ^= state[2];
// //         s3 ^= state[3];
// //       }
// //       operator()();
// //     }
// //
// //     state[0] = s0;
// //   state[1] = s1;
// //   state[2] = s2;
// //   state[3] = s3;
// // }
// //
// // /* This is the long-jump function for the generator. It is equivalent to
// //  2^192 calls to next(); it can be used to generate 2^64 starting points,
// //  from each of which jump() will generate 2^64 non-overlapping
// //  subsequences for parallel distributed computations. */
// // inline void xoshiro256plus::long_jump(void) {
// //   static const uint64_t LONG_JUMP[] = { 0x76e15d3efefdcbbf, 0xc5004e441c522fb3, 0x77710069854ee241, 0x39109bb02acbe635 };
// //
// //   uint64_t s0 = 0;
// //   uint64_t s1 = 0;
// //   uint64_t s2 = 0;
// //   uint64_t s3 = 0;
// //   for(unsigned int i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
// //     for(unsigned int b = 0; b < 64; b++) {
// //       if (LONG_JUMP[i] & uint64_t(1) << b) {
// //         s0 ^= state[0];
// //         s1 ^= state[1];
// //         s2 ^= state[2];
// //         s3 ^= state[3];
// //       }
// //       operator()();
// //     }
// //
// //     state[0] = s0;
// //   state[1] = s1;
// //   state[2] = s2;
// //   state[3] = s3;
// // }

#endif
