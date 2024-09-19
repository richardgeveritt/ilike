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
class MultiIndexParameters;

typedef boost::unordered_map< std::string, arma::mat>::iterator vector_parameter_iterator;
typedef boost::unordered_map< std::string, arma::mat>::const_iterator vector_parameter_const_iterator;
typedef boost::unordered_map< std::string, boost::spirit::hold_any>::iterator any_parameter_iterator;
typedef boost::unordered_map< std::string, boost::spirit::hold_any>::const_iterator any_parameter_const_iterator;

RCPP_EXPOSED_CLASS(MultiIndexParameters)

class MultiIndexParameters
{
  
public:
  
  MultiIndexParameters();
  
  virtual ~MultiIndexParameters();
  
  arma::mat& operator[](const std::string &variable);
  arma::mat operator[](const std::string &variable) const;
  arma::mat operator[](const std::vector<std::string> &variables) const;
  
  boost::spirit::hold_any& operator()(const std::string &variable);
  boost::spirit::hold_any operator()(const std::string &variable) const;
  
  MultiIndexParameters(const MultiIndexParameters &another);
  void operator=(const MultiIndexParameters &another);
  MultiIndexParameters* duplicate() const;
  
  //bool operator==(const MultiIndexParameters &another) const;
  //bool operator!=(const MultiIndexParameters &another) const;
  
  bool is_empty() const;
  
  arma::colvec get_vector() const;
  arma::colvec get_vector(const std::string &variable) const;
  arma::colvec get_vector(const std::vector<std::string> &variables) const;
  
  arma::rowvec get_row_vector() const;
  
  MultiIndexParameters merge(const MultiIndexParameters &another) const;
  
  vector_parameter_iterator vector_begin();
  vector_parameter_iterator vector_end();
  
  any_parameter_iterator any_begin();
  any_parameter_iterator any_end();
  
  vector_parameter_const_iterator vector_begin() const;
  vector_parameter_const_iterator vector_end() const;
  
  any_parameter_const_iterator any_begin() const;
  any_parameter_const_iterator any_end() const;
  
  size_t vector_size() const;
  size_t any_size() const;
  
  friend std::ostream& operator<<(std::ostream& os, const MultiIndexParameters &p);
  
protected:
  
  boost::unordered_map< std::string, arma::mat> vector_parameters;
  
  boost::unordered_map< std::string, boost::spirit::hold_any> any_parameters;
  
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
    return TRUE;
  else
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
  return concatenated_vector;
}

inline arma::colvec MultiIndexParameters::get_vector(const std::string &variable) const
{
  //arma::mat output = (*this)[variable];
  return arma::vectorise((*this)[variable]);
}

inline arma::colvec MultiIndexParameters::get_vector(const std::vector<std::string> &variables) const
{
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
  return concatenated_vector;
}

inline arma::rowvec MultiIndexParameters::get_row_vector() const
{
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
  
  return os;
}
}

#endif
