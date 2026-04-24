#ifndef MULTIINDEXPARAMETERS_H
#define MULTIINDEXPARAMETERS_H

#include <RcppArmadillo.h>
#include <RcppCommon.h>

#include <iostream>
#include <vector>
#include <string>
#include <boost/unordered_map.hpp>
#include <boost/spirit/home/support/detail/hold_any.hpp>


// At some point, we might need a more flexible container for data and/or parameters.
// At this point, make another class that includes this one.

namespace ilike
{
  /**
   * @file multiindex_parameters.h
   * @brief Defines the that class.
   *
   * Provides that functionality.
   *
   * @namespace ilike
   * @class that
   * @brief The that class.
   */


class MultiIndexParameters;

typedef boost::unordered_map< std::string, arma::mat>::iterator vector_parameter_iterator;
typedef boost::unordered_map< std::string, arma::mat>::const_iterator vector_parameter_const_iterator;
typedef boost::unordered_map< std::string, boost::spirit::hold_any>::iterator any_parameter_iterator;
typedef boost::unordered_map< std::string, boost::spirit::hold_any>::const_iterator any_parameter_const_iterator;

RCPP_EXPOSED_CLASS(MultiIndexParameters)

class MultiIndexParameters
{
  
public:
  
  /**
   * @brief Performs the multiindexparameters operation.
   */
  MultiIndexParameters();
  
  /**
   * @brief Performs the ~multiindexparameters operation.
   */
  virtual ~MultiIndexParameters();
  
  /**
   * @brief Performs the operator[] operation.
   *
   * @param variable The variable.
   *
   * @return The result.
   */
  arma::mat& operator[](const std::string &variable);
  /**
   * @brief Performs the operator[] operation.
   *
   * @param variable The variable.
   *
   * @return The result.
   */
  arma::mat operator[](const std::string &variable) const;
  /**
   * @brief Performs the operator[] operation.
   *
   * @param variables The variables.
   *
   * @return The result.
   */
  arma::mat operator[](const std::vector<std::string> &variables) const;
  
  /**
   * @brief Performs the operator() operation.
   *
   * @param variable The variable.
   *
   * @return The result.
   */
  boost::spirit::hold_any& operator()(const std::string &variable);
  /**
   * @brief Performs the operator() operation.
   *
   * @param variable The variable.
   *
   * @return The result.
   */
  boost::spirit::hold_any operator()(const std::string &variable) const;
  
  /**
   * @brief Performs the multiindexparameters operation.
   *
   * @param another The that instance to copy from.
   */
  MultiIndexParameters(const MultiIndexParameters &another);
  /**
   * @brief Assignment operator for that.
   *
   * @param another The that instance to copy from.
   */
  void operator=(const MultiIndexParameters &another);
  /**
   * @brief Creates a deep copy of this that object.
   *
   * @return The result.
   */
  MultiIndexParameters* duplicate() const;
  
  //bool operator==(const MultiIndexParameters &another) const;
  //bool operator!=(const MultiIndexParameters &another) const;
  
  /**
   * @brief Performs the is empty operation.
   *
   * @return The result.
   */
  bool is_empty() const;
  
  /**
   * @brief Returns the vector.
   *
   * @return The result.
   */
  arma::colvec get_vector() const;
  /**
   * @brief Returns the vector.
   *
   * @param variable The variable.
   *
   * @return The result.
   */
  arma::colvec get_vector(const std::string &variable) const;
  /**
   * @brief Returns the vector.
   *
   * @param variables The variables.
   *
   * @return The result.
   */
  arma::colvec get_vector(const std::vector<std::string> &variables) const;
  
  /**
   * @brief Returns the row vector.
   *
   * @return The result.
   */
  arma::rowvec get_row_vector() const;
  
  /**
   * @brief Performs the merge operation.
   *
   * @param another The that instance to copy from.
   *
   * @return The result.
   */
  MultiIndexParameters merge(const MultiIndexParameters &another) const;
  
  /**
   * @brief Performs the vector begin operation.
   *
   * @return The result.
   */
  vector_parameter_iterator vector_begin();
  /**
   * @brief Performs the vector end operation.
   *
   * @return The result.
   */
  vector_parameter_iterator vector_end();
  
  /**
   * @brief Performs the any begin operation.
   *
   * @return The result.
   */
  any_parameter_iterator any_begin();
  /**
   * @brief Performs the any end operation.
   *
   * @return The result.
   */
  any_parameter_iterator any_end();
  
  /**
   * @brief Performs the vector begin operation.
   *
   * @return The result.
   */
  vector_parameter_const_iterator vector_begin() const;
  /**
   * @brief Performs the vector end operation.
   *
   * @return The result.
   */
  vector_parameter_const_iterator vector_end() const;
  
  /**
   * @brief Performs the any begin operation.
   *
   * @return The result.
   */
  any_parameter_const_iterator any_begin() const;
  /**
   * @brief Performs the any end operation.
   *
   * @return The result.
   */
  any_parameter_const_iterator any_end() const;
  
  /**
   * @brief Performs the vector size operation.
   *
   * @return The result.
   */
  size_t vector_size() const;
  /**
   * @brief Performs the any size operation.
   *
   * @return The result.
   */
  size_t any_size() const;
  
  /**
   * @brief Performs the operator<< operation.
   *
   * @param os The os.
   * @param p The p.
   *
   * @return The result.
   */
  friend std::ostream& operator<<(std::ostream& os, const MultiIndexParameters &p);
  
protected:
  
  boost::unordered_map< std::string, arma::mat> vector_parameters;
  
  boost::unordered_map< std::string, boost::spirit::hold_any> any_parameters;
  
  /**
   * @brief Copies the state of another that into this object.
   *
   * @param another The that instance to copy from.
   */
  void make_copy(const MultiIndexParameters &another);
  
};

//typedef MultiIndexParameters Data;

//RCPP_EXPOSED_WRAP(MultiIndexParameters);

#include <RcppArmadillo.h>

inline MultiIndexParameters::MultiIndexParameters()
{
  this->vector_parameters.clear();
  this->any_parameters.clear();
}

inline MultiIndexParameters::~MultiIndexParameters()
{
  
}

//Copy constructor for the MultiIndexParameters class.
inline MultiIndexParameters::MultiIndexParameters(const MultiIndexParameters &another)
{
  this->make_copy(another);
}

inline void MultiIndexParameters::operator=(const MultiIndexParameters &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  this->make_copy(another);
}

inline MultiIndexParameters* MultiIndexParameters::duplicate() const
{
  return( new MultiIndexParameters(*this));
}

inline bool MultiIndexParameters::is_empty() const
{
  if ( (this->vector_parameters.size()==0) && (this->any_parameters.size()==0) )
    /** @brief The true. */
    return TRUE;
  else
    /** @brief The false. */
    return FALSE;
}

/*
 inline bool MultiIndexParameters::operator==(const MultiIndexParameters &another) const
 {
 if ((*this)!=another)
 {
 return false;
 }
 else
 {
 return true;
 }
 }
 
 inline bool MultiIndexParameters::operator!=(const MultiIndexParameters &another) const
 {
 if (this->vector_parameters != another.vector_parameters)
 {
 return true;
 }
 
 if (this->any_parameters != another.any_parameters)
 {
 return true;
 }
 
 return false;
 }
 */

inline arma::colvec MultiIndexParameters::get_vector() const
{
  /** @brief The concatenated vector. */
  arma::colvec concatenated_vector;
  for (vector_parameter_const_iterator i=this->vector_begin();
       i!=this->vector_end();
       ++i)
  {
    if (i==this->vector_begin())
      concatenated_vector = arma::vectorise(this->vector_begin()->second);
    else
      concatenated_vector = join_cols(concatenated_vector,arma::vectorise(i->second));
  }
  /** @brief The concatenated vector. */
  return concatenated_vector;
}

inline arma::colvec MultiIndexParameters::get_vector(const std::string &variable) const
{
  //arma::mat output = (*this)[variable];
  return arma::vectorise((*this)[variable]);
}

inline arma::colvec MultiIndexParameters::get_vector(const std::vector<std::string> &variables) const
{
  /** @brief The concatenated vector. */
  arma::colvec concatenated_vector;
  for (auto i=variables.begin();
       i!=variables.end();
       ++i)
  {
    if (i==variables.begin())
      concatenated_vector = this->get_vector(*i);
    else
      concatenated_vector = join_cols(concatenated_vector,this->get_vector(*i));
  }
  /** @brief The concatenated vector. */
  return concatenated_vector;
}

inline arma::rowvec MultiIndexParameters::get_row_vector() const
{
  /** @brief The concatenated vector. */
  arma::rowvec concatenated_vector;
  for (vector_parameter_const_iterator i=this->vector_begin();
       i!=this->vector_end();
       ++i)
  {
    if (i==this->vector_begin())
      concatenated_vector = arma::vectorise(this->vector_begin()->second);
    else
      concatenated_vector = join_rows(concatenated_vector,arma::vectorise(i->second));
  }
  /** @brief The concatenated vector. */
  return concatenated_vector;
}

inline void MultiIndexParameters::make_copy(const MultiIndexParameters &another)
{
  this->vector_parameters = another.vector_parameters;
  this->any_parameters = another.any_parameters;
}

inline arma::mat& MultiIndexParameters::operator[](const std::string &variable)
{
  return(this->vector_parameters[variable]);
}

inline arma::mat MultiIndexParameters::operator[](const std::string &variable) const
{
  boost::unordered_map< std::string, arma::mat>::const_iterator found =  this->vector_parameters.find(variable);
  
  if (found != this->vector_parameters.end())
    return(found->second);
  else
    return arma::mat();
}

inline arma::mat MultiIndexParameters::operator[](const std::vector<std::string> &variables) const
{
  /** @brief The concatenated matrix. */
  arma::mat concatenated_matrix;
  for (std::vector<std::string>::const_iterator i=variables.begin();
       i!=variables.end();
       ++i)
  {
    boost::unordered_map< std::string, arma::mat>::const_iterator found =  this->vector_parameters.find(*i);
    if (i==variables.begin())
    {
      if (found != this->vector_parameters.end())
        concatenated_matrix = found->second;
      else
        Rcpp::stop("MultiIndexParameters::operator[]: variable not found in MultiIndexParameters.");
    }
    else
    {
      if (found != this->vector_parameters.end())
        concatenated_matrix = join_cols(concatenated_matrix,found->second);
      else
        Rcpp::stop("MultiIndexParameters::operator[]: variable not found in MultiIndexParameters.");
    }
  }
  
  /** @brief The concatenated matrix. */
  return concatenated_matrix;
}

inline boost::spirit::hold_any& MultiIndexParameters::operator()(const std::string &variable)
{
  return(this->any_parameters[variable]);
}

inline boost::spirit::hold_any MultiIndexParameters::operator()(const std::string &variable) const
{
  boost::unordered_map< std::string, boost::spirit::hold_any>::const_iterator found = this->any_parameters.find(variable);
  
  if (found != this->any_parameters.end())
    return(found->second);
  else
    Rcpp::stop("MultiIndexParameters::operator(): variable not found in MultiIndexParameters.");
}

inline MultiIndexParameters MultiIndexParameters::merge(const MultiIndexParameters &another) const
{
  MultiIndexParameters output = *this;
  
  output.vector_parameters.insert(another.vector_parameters.begin(), another.vector_parameters.end());
  
  output.any_parameters.insert(another.any_parameters.begin(), another.any_parameters.end());
  
  /** @brief The output. */
  return output;
}

inline vector_parameter_iterator MultiIndexParameters::vector_begin()
{
  return this->vector_parameters.begin();
}

inline vector_parameter_iterator MultiIndexParameters::vector_end()
{
  return this->vector_parameters.end();
}

inline any_parameter_iterator MultiIndexParameters::any_begin()
{
  return this->any_parameters.begin();
}

inline any_parameter_iterator MultiIndexParameters::any_end()
{
  return this->any_parameters.end();
}

inline vector_parameter_const_iterator MultiIndexParameters::vector_begin() const
{
  return this->vector_parameters.begin();
}

inline vector_parameter_const_iterator MultiIndexParameters::vector_end() const
{
  return this->vector_parameters.end();
}

inline any_parameter_const_iterator MultiIndexParameters::any_begin() const
{
  return this->any_parameters.begin();
}

inline any_parameter_const_iterator MultiIndexParameters::any_end() const
{
  return this->any_parameters.end();
}

inline size_t MultiIndexParameters::vector_size() const
{
  return this->vector_parameters.size();
}

inline size_t MultiIndexParameters::any_size() const
{
  return this->any_parameters.size();
}

inline std::ostream& operator<<(std::ostream& os, const MultiIndexParameters &p)
{
  boost::unordered_map< std::string, arma::mat>::const_iterator it;
  
  for (it=p.vector_parameters.begin();it!=p.vector_parameters.end();++it)
  {
    os << it->first << ":";
    for (arma::mat::const_iterator i=it->second.begin(); i!=it->second.end(); ++i)
    {
      if (i==it->second.begin())
        os << *i;
      else
        os << "," << *i;
    }
    os << std::endl;
  }
  
  boost::unordered_map< std::string, boost::spirit::hold_any>::const_iterator it2;
  
  for (it2=p.any_parameters.begin();it2!=p.any_parameters.end();++it2)
  {
    os << it2->first << ":";
    os << it2->second;
    os << std::endl;
  }
  
  /** @brief The os. */
  return os;
}
}

#endif
