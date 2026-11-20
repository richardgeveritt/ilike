// Automatic differentiation for log densities written by the user as
// templates on the number type (gradient = autodiff, dev/plan.md 3.10
// stages 2 and 4). Included only in models whose C++ code uses ParametersT,
// so the package's own sources do not depend on it.
//
//   template <class T>
//   T log_lik(const ParametersT<T> &parameters, const Data &data) { ... }
//
// ilike instantiates the function with double, through a ParametersT<double>
// that is a view on the Parameters (the log density then costs what it costs
// when written with Parameters), and for the gradient with one of:
// - Var (reverse mode, the default): one evaluation records each operation
//   on a tape, and one sweep back along the tape gives the derivatives with
//   respect to every element of the variable;
// - Dual (forward mode, autodiff(forward)): a number carrying a value and its
//   derivative in one direction, one evaluation per element.

#ifndef ILIKE_AUTODIFF_H
#define ILIKE_AUTODIFF_H

#include <cmath>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>
#include <boost/math/special_functions/digamma.hpp>
#include "parameters.h"

namespace ilike {

// --- forward mode: dual numbers ---------------------------------------------------

// a value and its derivative in one direction: arithmetic on Duals applies
// the chain rule to the derivative alongside the value
struct Dual
{
  double value;
  double derivative;
  Dual(double value_in = 0.0, double derivative_in = 0.0) : value(value_in), derivative(derivative_in) {}
  Dual &operator+=(const Dual &b) { this->value += b.value; this->derivative += b.derivative; return *this; }
  Dual &operator-=(const Dual &b) { this->value -= b.value; this->derivative -= b.derivative; return *this; }
  Dual &operator*=(const Dual &b) { this->derivative = this->derivative * b.value + this->value * b.derivative; this->value *= b.value; return *this; }
  Dual &operator/=(const Dual &b) { this->derivative = (this->derivative * b.value - this->value * b.derivative) / (b.value * b.value); this->value /= b.value; return *this; }
};

inline Dual operator+(const Dual &a) { return a; }
inline Dual operator-(const Dual &a) { return Dual(-a.value, -a.derivative); }
inline Dual operator+(const Dual &a, const Dual &b) { return Dual(a.value + b.value, a.derivative + b.derivative); }
inline Dual operator+(const Dual &a, double b) { return Dual(a.value + b, a.derivative); }
inline Dual operator+(double a, const Dual &b) { return Dual(a + b.value, b.derivative); }
inline Dual operator-(const Dual &a, const Dual &b) { return Dual(a.value - b.value, a.derivative - b.derivative); }
inline Dual operator-(const Dual &a, double b) { return Dual(a.value - b, a.derivative); }
inline Dual operator-(double a, const Dual &b) { return Dual(a - b.value, -b.derivative); }
inline Dual operator*(const Dual &a, const Dual &b) { return Dual(a.value * b.value, a.derivative * b.value + a.value * b.derivative); }
inline Dual operator*(const Dual &a, double b) { return Dual(a.value * b, a.derivative * b); }
inline Dual operator*(double a, const Dual &b) { return Dual(a * b.value, a * b.derivative); }
inline Dual operator/(const Dual &a, const Dual &b) { return Dual(a.value / b.value, (a.derivative * b.value - a.value * b.derivative) / (b.value * b.value)); }
inline Dual operator/(const Dual &a, double b) { return Dual(a.value / b, a.derivative / b); }
inline Dual operator/(double a, const Dual &b) { return Dual(a / b.value, -a * b.derivative / (b.value * b.value)); }

// functions: f(a) has derivative f'(a.value) a.derivative
inline Dual dual_function(const Dual &a, double value, double slope) { return Dual(value, slope * a.derivative); }

// --- reverse mode: numbers recorded on a tape -------------------------------------

// one operation: the result's parents (-1: none) and its partial derivative
// with respect to each
struct TapeNode
{
  int a;
  int b;
  double da;
  double db;
};

// the tape of the current thread (its capacity is kept between gradients)
inline std::vector<TapeNode> &reverse_tape()
{
  thread_local std::vector<TapeNode> tape;
  return tape;
}

// a value, its place on the tape (-1 for a constant, which is not recorded)
// and the tape (carried, so that an operation does not look up the thread's
// tape, a call on some platforms, for every operation)
struct Var
{
  double value;
  int index;
  std::vector<TapeNode> *tape;
  Var(double value_in = 0.0) : value(value_in), index(-1), tape(nullptr) {}
  Var(double value_in, int index_in, std::vector<TapeNode> *tape_in) : value(value_in), index(index_in), tape(tape_in) {}
  Var &operator+=(const Var &b);
  Var &operator-=(const Var &b);
  Var &operator*=(const Var &b);
  Var &operator/=(const Var &b);
};

inline Var var_unary(double value, const Var &a, double da)
{
  if (a.index < 0)
    return Var(value);
  a.tape->push_back(TapeNode{a.index, -1, da, 0.0});
  return Var(value, (int)a.tape->size() - 1, a.tape);
}

inline Var var_binary(double value, const Var &a, double da, const Var &b, double db)
{
  if (a.index < 0)
    return var_unary(value, b, db);
  if (b.index < 0)
    return var_unary(value, a, da);
  a.tape->push_back(TapeNode{a.index, b.index, da, db});
  return Var(value, (int)a.tape->size() - 1, a.tape);
}

inline Var operator+(const Var &a) { return a; }
inline Var operator-(const Var &a) { return var_unary(-a.value, a, -1.0); }
inline Var operator+(const Var &a, const Var &b) { return var_binary(a.value + b.value, a, 1.0, b, 1.0); }
inline Var operator+(const Var &a, double b) { return var_unary(a.value + b, a, 1.0); }
inline Var operator+(double a, const Var &b) { return var_unary(a + b.value, b, 1.0); }
inline Var operator-(const Var &a, const Var &b) { return var_binary(a.value - b.value, a, 1.0, b, -1.0); }
inline Var operator-(const Var &a, double b) { return var_unary(a.value - b, a, 1.0); }
inline Var operator-(double a, const Var &b) { return var_unary(a - b.value, b, -1.0); }
inline Var operator*(const Var &a, const Var &b) { return var_binary(a.value * b.value, a, b.value, b, a.value); }
inline Var operator*(const Var &a, double b) { return var_unary(a.value * b, a, b); }
inline Var operator*(double a, const Var &b) { return var_unary(a * b.value, b, a); }
inline Var operator/(const Var &a, const Var &b) { return var_binary(a.value / b.value, a, 1.0 / b.value, b, -a.value / (b.value * b.value)); }
inline Var operator/(const Var &a, double b) { return var_unary(a.value / b, a, 1.0 / b); }
inline Var operator/(double a, const Var &b) { return var_unary(a / b.value, b, -a / (b.value * b.value)); }
inline Var &Var::operator+=(const Var &b) { *this = *this + b; return *this; }
inline Var &Var::operator-=(const Var &b) { *this = *this - b; return *this; }
inline Var &Var::operator*=(const Var &b) { *this = *this * b; return *this; }
inline Var &Var::operator/=(const Var &b) { *this = *this / b; return *this; }

inline Var dual_function(const Var &a, double value, double slope) { return var_unary(value, a, slope); }

// --- what the two types share ------------------------------------------------------

inline double value_of(double a) { return a; }
inline double value_of(const Dual &a) { return a.value; }
inline double value_of(const Var &a) { return a.value; }

// comparisons compare values (a branch on a value is constant near it)
#define ILIKE_AD_COMPARISON(N, op) \
  inline bool operator op(const N &a, const N &b) { return a.value op b.value; } \
  inline bool operator op(const N &a, double b) { return a.value op b; } \
  inline bool operator op(double a, const N &b) { return a op b.value; }
#define ILIKE_AD_COMPARISONS(N) \
  ILIKE_AD_COMPARISON(N, <) ILIKE_AD_COMPARISON(N, >) ILIKE_AD_COMPARISON(N, <=) \
  ILIKE_AD_COMPARISON(N, >=) ILIKE_AD_COMPARISON(N, ==) ILIKE_AD_COMPARISON(N, !=)
ILIKE_AD_COMPARISONS(Dual)
ILIKE_AD_COMPARISONS(Var)
#undef ILIKE_AD_COMPARISONS
#undef ILIKE_AD_COMPARISON

// functions: f(a) = f(a.value), with slope f'(a.value)
#define ILIKE_AD_FUNCTIONS(N) \
  inline N exp(const N &a) { double e = std::exp(a.value); return dual_function(a, e, e); } \
  inline N log(const N &a) { return dual_function(a, std::log(a.value), 1.0 / a.value); } \
  inline N sqrt(const N &a) { double s = std::sqrt(a.value); return dual_function(a, s, 0.5 / s); } \
  inline N log1p(const N &a) { return dual_function(a, std::log1p(a.value), 1.0 / (1.0 + a.value)); } \
  inline N expm1(const N &a) { return dual_function(a, std::expm1(a.value), std::exp(a.value)); } \
  inline N sin(const N &a) { return dual_function(a, std::sin(a.value), std::cos(a.value)); } \
  inline N cos(const N &a) { return dual_function(a, std::cos(a.value), -std::sin(a.value)); } \
  inline N tan(const N &a) { double t = std::tan(a.value); return dual_function(a, t, 1.0 + t * t); } \
  inline N sinh(const N &a) { return dual_function(a, std::sinh(a.value), std::cosh(a.value)); } \
  inline N cosh(const N &a) { return dual_function(a, std::cosh(a.value), std::sinh(a.value)); } \
  inline N tanh(const N &a) { double t = std::tanh(a.value); return dual_function(a, t, 1.0 - t * t); } \
  inline N atan(const N &a) { return dual_function(a, std::atan(a.value), 1.0 / (1.0 + a.value * a.value)); } \
  inline N abs(const N &a) { return (a.value < 0.0) ? -a : a; } /* (one-sided at 0) */ \
  inline N fabs(const N &a) { return abs(a); } \
  inline N lgamma(const N &a) { return dual_function(a, std::lgamma(a.value), boost::math::digamma(a.value)); } \
  inline N pow(const N &a, double b) \
  { \
    if (b == 0.0) \
      return N(1.0); \
    return dual_function(a, std::pow(a.value, b), b * std::pow(a.value, b - 1.0)); \
  } \
  inline N pow(double a, const N &b) { double p = std::pow(a, b.value); return dual_function(b, p, p * std::log(a)); }
ILIKE_AD_FUNCTIONS(Dual)
ILIKE_AD_FUNCTIONS(Var)
#undef ILIKE_AD_FUNCTIONS

// pow with both arguments varying (the log term skipped where the exponent
// does not vary, so a negative base with a constant exponent is fine)
inline Dual pow(const Dual &a, const Dual &b)
{
  double p = std::pow(a.value, b.value);
  double from_a = (a.derivative == 0.0) ? 0.0 : b.value * std::pow(a.value, b.value - 1.0) * a.derivative;
  double from_b = (b.derivative == 0.0) ? 0.0 : p * std::log(a.value) * b.derivative;
  return Dual(p, from_a + from_b);
}
inline Var pow(const Var &a, const Var &b)
{
  if (b.index < 0)
    return pow(a, b.value);
  if (a.index < 0)
    return pow(a.value, b);
  double p = std::pow(a.value, b.value);
  return var_binary(p, a, b.value * std::pow(a.value, b.value - 1.0), b, p * std::log(a.value));
}

// --- log densities (as the double versions in distributions.h) ---------------------

// the number type of a call whose arguments are doubles or one autodiff type
// (no type, so the double versions are chosen, when no argument is one)
template <class T> struct is_autodiff : std::false_type {};
template <> struct is_autodiff<Dual> : std::true_type {};
template <> struct is_autodiff<Var> : std::true_type {};
struct autodiff_none {}; // (an argument that is neither a number nor the autodiff type)
template <class N, class T, class = void> struct autodiff_join { typedef autodiff_none type; };
template <class T> struct autodiff_join<void, T, typename std::enable_if<std::is_arithmetic<T>::value>::type> { typedef void type; };
template <class T> struct autodiff_join<void, T, typename std::enable_if<is_autodiff<T>::value>::type> { typedef T type; };
template <class N, class T> struct autodiff_join<N, T, typename std::enable_if<is_autodiff<N>::value && std::is_arithmetic<T>::value>::type> { typedef N type; };
template <class N> struct autodiff_join<N, N, typename std::enable_if<is_autodiff<N>::value>::type> { typedef N type; };
template <class N, class... T> struct autodiff_fold { typedef N type; };
template <class N, class T, class... R> struct autodiff_fold<N, T, R...>
  : autodiff_fold<typename autodiff_join<N, typename std::decay<T>::type>::type, R...> {};
// (type only when the arguments make an autodiff type: otherwise the
// template is set aside, without an error, and the double version is used)
template <class F, class = void> struct autodiff_found {};
template <class F> struct autodiff_found<F, typename std::enable_if<is_autodiff<F>::value>::type> { typedef F type; };
template <class... T> struct autodiff_number : autodiff_found<typename autodiff_fold<void, T...>::type> {};
template <class... T> using autodiff_t = typename autodiff_number<T...>::type;

template <class A, class N = autodiff_t<A>>
N dnorm(const A &x_in)
{
  N x(x_in);
  return -0.5 * std::log(2.0 * M_PI) - 0.5 * x * x;
}
template <class A, class B, class C, class N = autodiff_t<A, B, C>>
N dnorm(const A &x_in, const B &mean_in, const C &sd_in)
{
  N x(x_in), mean(mean_in), sd(sd_in);
  if (sd.value < 0.0)
    return N(NAN);
  if (sd.value == 0.0)
    return N((x.value == mean.value) ? arma::datum::inf : -arma::datum::inf);
  N z = (x - mean) / sd;
  return -log(sd) - 0.5 * std::log(2.0 * M_PI) - 0.5 * z * z;
}
template <class A, class B, class C, class N = autodiff_t<A, B, C>>
N dlnorm(const A &x_in, const B &meanlog_in, const C &sdlog_in)
{
  N x(x_in), meanlog(meanlog_in), sdlog(sdlog_in);
  if (sdlog.value < 0.0)
    return N(NAN);
  if (sdlog.value == 0.0)
    return N((x.value == meanlog.value) ? arma::datum::inf : -arma::datum::inf);
  N z = (log(x) - meanlog) / sdlog;
  return -log(x) - log(sdlog) - 0.5 * std::log(2.0 * M_PI) - 0.5 * z * z;
}
template <class A, class B, class N = autodiff_t<A, B>>
N dexp(const A &x_in, const B &rate_in)
{
  N x(x_in), rate(rate_in);
  if (rate.value <= 0.0)
    return N(NAN);
  if (x.value < 0.0)
    return N(-arma::datum::inf);
  return log(rate) - rate * x;
}
template <class A, class B, class C, class N = autodiff_t<A, B, C>>
N dgamma(const A &x_in, const B &shape_in, const C &rate_in)
{
  N x(x_in), shape(shape_in), rate(rate_in);
  if ((shape.value <= 0.0) || (rate.value <= 0.0))
    return N(NAN);
  if (x.value < 0.0)
    return N(-arma::datum::inf);
  return -lgamma(shape) + shape * log(rate) + (shape - 1.0) * log(x) - rate * x;
}
template <class A, class B, class N = autodiff_t<A, B>>
N dpois(const A &x_in, const B &rate_in)
{
  double x = value_of(x_in); // (a count: not differentiated)
  N rate(rate_in);
  if (rate.value <= 0.0)
    return N(NAN);
  if (x < 0.0)
    return N(-arma::datum::inf);
  return x * log(rate) - rate - std::lgamma(x + 1.0);
}
template <class A, class N = autodiff_t<A>>
N dunif(const A &x_in)
{
  double x = value_of(x_in);
  return N(((x < 0.0) || (x > 1.0)) ? -arma::datum::inf : 0.0);
}
template <class A, class B, class C, class N = autodiff_t<A, B, C>>
N dunif(const A &x_in, const B &lower_in, const C &upper_in)
{
  N lower(lower_in), upper(upper_in);
  double x = value_of(x_in);
  if ((x < lower.value) || (x > upper.value))
    return N(-arma::datum::inf);
  return -log(upper - lower);
}

// --- parameters seen by a templated log density -----------------------------------

template <class T> class ParametersT;

// with double: a view on the Parameters, giving their matrices as they are
template <> class ParametersT<double>
{
public:
  explicit ParametersT(const Parameters &parameters_in) : parameters(parameters_in) {}
  const arma::mat &operator[](const std::string &variable) const { return this->parameters[variable]; }
private:
  const Parameters &parameters;
};

// the elements of a variable as Duals: the seeded element has derivative 1,
// the others 0 (the elements of the other variables too)
class DualMatrix
{
public:
  DualMatrix(const arma::mat &values_in, arma::uword seeded_in)
    : n_rows(values_in.n_rows), n_cols(values_in.n_cols), n_elem(values_in.n_elem), values(&values_in), seeded(seeded_in) {}
  Dual operator[](arma::uword i) const { return Dual((*this->values)(i), (i == this->seeded) ? 1.0 : 0.0); }
  Dual operator()(arma::uword i) const { return (*this)[i]; }
  Dual operator()(arma::uword row, arma::uword column) const { return (*this)[row + column * this->n_rows]; }
  const arma::uword n_rows;
  const arma::uword n_cols;
  const arma::uword n_elem;
private:
  const arma::mat *values;
  arma::uword seeded;
};

// with Dual: the derivative with respect to one element of one variable
template <> class ParametersT<Dual>
{
public:
  ParametersT(const Parameters &parameters_in, const std::string &variable_in, arma::uword element_in)
    : parameters(parameters_in), variable(variable_in), element(element_in) {}
  DualMatrix operator[](const std::string &name) const
  {
    return DualMatrix(this->parameters[name], (name == this->variable) ? this->element : std::numeric_limits<arma::uword>::max());
  }
private:
  const Parameters &parameters;
  const std::string &variable;
  arma::uword element;
};

// the elements of a variable as Vars: those of the variable differentiated
// are the tape's leaves 0, 1, ...; the others are constants
class VarMatrix
{
public:
  VarMatrix(const arma::mat &values_in, std::vector<TapeNode> *tape_in)
    : n_rows(values_in.n_rows), n_cols(values_in.n_cols), n_elem(values_in.n_elem), values(&values_in), tape(tape_in) {}
  Var operator[](arma::uword i) const { return (this->tape != nullptr) ? Var((*this->values)(i), (int)i, this->tape) : Var((*this->values)(i)); }
  Var operator()(arma::uword i) const { return (*this)[i]; }
  Var operator()(arma::uword row, arma::uword column) const { return (*this)[row + column * this->n_rows]; }
  const arma::uword n_rows;
  const arma::uword n_cols;
  const arma::uword n_elem;
private:
  const arma::mat *values;
  std::vector<TapeNode> *tape; // (nullptr: not the variable differentiated)
};

// with Var: the derivatives with respect to every element of one variable
template <> class ParametersT<Var>
{
public:
  ParametersT(const Parameters &parameters_in, const std::string &variable_in, std::vector<TapeNode> *tape_in)
    : parameters(parameters_in), variable(variable_in), tape(tape_in) {}
  VarMatrix operator[](const std::string &name) const
  {
    return VarMatrix(this->parameters[name], (name == this->variable) ? this->tape : nullptr);
  }
private:
  const Parameters &parameters;
  const std::string &variable;
  std::vector<TapeNode> *tape;
};

// --- gradients ---------------------------------------------------------------------

// forward mode: one Dual evaluation of f (called as f(parameters_T)) per element
template <class F>
arma::mat forward_gradient(F f, const std::string &variable, const Parameters &parameters)
{
  const arma::mat &x = parameters[variable];
  arma::mat gradient(x.n_rows, x.n_cols);
  for (arma::uword k = 0; k < x.n_elem; ++k)
    gradient(k) = f(ParametersT<Dual>(parameters, variable, k)).derivative;
  return gradient;
}

// reverse mode: one Var evaluation of f, recording its operations, then one
// sweep back along the tape, each node passing its adjoint (the derivative
// of the result with respect to it) to its parents through its partials
template <class F>
arma::mat reverse_gradient(F f, const std::string &variable, const Parameters &parameters)
{
  const arma::mat &x = parameters[variable];
  std::vector<TapeNode> &tape = reverse_tape();
  tape.clear();
  for (arma::uword k = 0; k < x.n_elem; ++k)
    tape.push_back(TapeNode{-1, -1, 0.0, 0.0}); // (the leaves)
  Var result = f(ParametersT<Var>(parameters, variable, &tape));
  arma::mat gradient(x.n_rows, x.n_cols, arma::fill::zeros);
  if (result.index < 0)
    return gradient; // (the result does not depend on the variable)
  thread_local std::vector<double> adjoint;
  adjoint.assign(tape.size(), 0.0);
  adjoint[result.index] = 1.0;
  for (int i = result.index; i >= (int)x.n_elem; --i)
  {
    const double a = adjoint[i];
    if (a == 0.0)
      continue;
    const TapeNode &node = tape[i];
    if (node.a >= 0)
      adjoint[node.a] += node.da * a;
    if (node.b >= 0)
      adjoint[node.b] += node.db * a;
  }
  for (arma::uword k = 0; k < x.n_elem; ++k)
    gradient(k) = adjoint[k];
  return gradient;
}

}

#endif
