// The R binding's side of the recipe (dev/plan.md 3.14, I1): the R lists made
// by compile() and passed to the do_* functions are converted here, once per
// run, into RecipeValue trees, which are all the core reads.

#include <RcppArmadillo.h>
#include <climits>
#include "recipe_value.h"
#include "recipe_from_r.h"
#include "builtin_registry.h"
#include "ilike_language.h"

using namespace Rcpp;

namespace ilike
{

static RecipeValue convert(SEXP x, const std::string &name)
{
  switch (TYPEOF(x))
  {
    case NILSXP:
      return RecipeValue();
    case VECSXP:
    {
      RecipeValue output = RecipeValue::list();
      SEXP names = Rf_getAttrib(x, R_NamesSymbol);
      R_xlen_t n = Rf_xlength(x);
      for (R_xlen_t i = 0; i < n; ++i)
      {
        std::string key = "";
        if (!Rf_isNull(names))
        {
          SEXP k = STRING_ELT(names, i);
          if (k != NA_STRING)
            key = CHAR(k);
        }
        output.push_back(convert(VECTOR_ELT(x, i), key), key);
      }
      return output;
    }
    case REALSXP:
    {
      if (Rf_isMatrix(x))
        return RecipeValue(as<arma::mat>(x));
      R_xlen_t n = Rf_xlength(x);
      std::vector<double> values(REAL(x), REAL(x) + n);
      return RecipeValue(values);
    }
    case INTSXP:
    {
      if (Rf_isMatrix(x))
        return RecipeValue(as<arma::mat>(x));
      if (Rf_isFactor(x))
        return RecipeValue::other("factor");
      R_xlen_t n = Rf_xlength(x);
      return RecipeValue::integers(std::vector<int>(INTEGER(x), INTEGER(x) + n));
    }
    case LGLSXP:
    {
      R_xlen_t n = Rf_xlength(x);
      return RecipeValue::logicals(std::vector<int>(LOGICAL(x), LOGICAL(x) + n));
    }
    case STRSXP:
    {
      R_xlen_t n = Rf_xlength(x);
      std::vector<std::string> values;
      values.reserve(n);
      for (R_xlen_t i = 0; i < n; ++i)
      {
        SEXP s = STRING_ELT(x, i);
        values.push_back(s == NA_STRING ? std::string("NA") : std::string(CHAR(s)));
      }
      return RecipeValue(values);
    }
    case EXTPTRSXP:
    {
      // a function compiled from the recipe: XPtr<FunctionPointer> holds the
      // address of the C++ function pointer
      FunctionReference reference;
      reference.language = "cpp";
      reference.name = name;
      reference.pointer = R_ExternalPtrAddr(x);
      return RecipeValue(reference);
    }
    case CLOSXP:
      return RecipeValue::other("R function");
    case ENVSXP:
      return RecipeValue::other("R environment");
    default:
      return RecipeValue::other("R object of type " + std::string(Rf_type2char(TYPEOF(x))));
  }
}

RecipeValue recipe_from_R(SEXP x)
{
  return convert(x, "");
}

SEXP recipe_to_R(const RecipeValue &x)
{
  switch (x.kind())
  {
    case RecipeValue::Kind::Null:
      return R_NilValue;
    case RecipeValue::Kind::Logical:
    {
      std::vector<int> values = x.as_ints();
      LogicalVector output(values.size());
      for (size_t i = 0; i < values.size(); ++i)
        output[i] = values[i];
      return output;
    }
    case RecipeValue::Kind::Integer:
      return wrap(x.as_ints());
    case RecipeValue::Kind::Real:
      return wrap(x.as_doubles());
    case RecipeValue::Kind::String:
      return wrap(x.as_strings());
    case RecipeValue::Kind::Matrix:
      return wrap(x.as_matrix());
    case RecipeValue::Kind::List:
    {
      List output(x.size());
      std::vector<std::string> keys = x.names();
      bool named = false;
      for (size_t i = 0; i < x.size(); ++i)
      {
        output[i] = recipe_to_R(x[i]);
        if ((i < keys.size()) && !keys[i].empty())
          named = true;
      }
      if (named)
        output.attr("names") = wrap(keys);
      return output;
    }
    case RecipeValue::Kind::Function:
      return wrap(x.as_function().language + "::" + x.as_function().name);
    case RecipeValue::Kind::Other:
      return R_NilValue;
  }
  return R_NilValue;
}

}

// [[Rcpp::export]]
std::string recipe_json(SEXP recipe)
{
  return ilike::recipe_from_R(recipe).to_json();
}

// The registry of built-ins (dev/plan.md 3.14, I2), for builtin_table().
// [[Rcpp::export]]
SEXP builtin_registry_table()
{
  return ilike::recipe_to_R(ilike::builtin_table());
}

// The syntax tree of a file in the ilike language (dev/ilike-language.md,
// I3a), for inspection and tests: as an R list, or as JSON.
// [[Rcpp::export]]
SEXP ilike_language_syntax_tree(const std::string &text, const std::string &file_name, bool json)
{
  ilike::RecipeValue tree = ilike::parse_ilike_language(text, file_name);
  if (json)
    return Rcpp::wrap(tree.to_json());
  return ilike::recipe_to_R(tree);
}

// [[Rcpp::export]]
bool is_ilike_language_text(const std::string &text)
{
  return ilike::is_ilike_language(text);
}

// Parse and lower the files of a recipe in the ilike language (I3b): the
// recipe with placeholders for the user's functions, the algorithms, the
// constants' defaults and the C++ code to compile (R/ilike_language.R).
// [[Rcpp::export]]
SEXP ilike_language_lower(const std::vector<std::string> &texts, const std::vector<std::string> &file_names,
                          SEXP constant_values = R_NilValue)
{
  ilike::RecipeValue trees = ilike::RecipeValue::list();
  for (size_t i = 0; i < texts.size(); ++i)
    trees.push_back(ilike::parse_ilike_language(texts[i], file_names[i]));
  ilike::RecipeValue values = Rf_isNull(constant_values) ? ilike::RecipeValue() : ilike::recipe_from_R(constant_values);
  return ilike::recipe_to_R(ilike::lower_ilike_language(trees, values));
}

// [[Rcpp::export]]
std::string ilike_find_cpp_function_signature(const std::string &code, const std::string &name)
{
  return ilike::find_cpp_function_signature(code, name);
}
