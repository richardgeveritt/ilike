#ifndef DATA_H
#define DATA_H

#include <RcppArmadilloForward.h>
#include <RcppCommon.h>

#include <iostream>
#include <vector>
#include <string>
#include <boost/unordered_map.hpp>
#include <boost/spirit/home/support/detail/hold_any.hpp>


// At some point, we might need a more flexible container for data and/or data.
// At this point, make another class that includes this one.

class Data;

typedef boost::unordered_map< std::string, arma::mat>::iterator vector_parameter_iterator;
typedef boost::unordered_map< std::string, arma::mat>::const_iterator vector_parameter_const_iterator;
typedef boost::unordered_map< std::string, boost::spirit::hold_any>::iterator any_parameter_iterator;
typedef boost::unordered_map< std::string, boost::spirit::hold_any>::const_iterator any_parameter_const_iterator;

RCPP_EXPOSED_CLASS(Data)

class Data
{

public:

  Data();

  virtual ~Data();

  arma::mat& operator[](const std::string &variable);
  arma::mat operator[](const std::string &variable) const;
  arma::mat operator[](const std::vector<std::string> &variables) const;
  
  boost::spirit::hold_any& operator()(const std::string &variable);
  boost::spirit::hold_any operator()(const std::string &variable) const;

  Data(const Data &another);
  void operator=(const Data &another);
  Data* duplicate() const;
  
  //bool operator==(const Data &another) const;
  //bool operator!=(const Data &another) const;
  
  bool is_empty() const;
  
  arma::colvec get_vector() const;
  arma::colvec get_vector(const std::string &variable) const;
  arma::colvec get_vector(const std::vector<std::string> &variables) const;
  
  Data merge(const Data &another) const;
  
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

  friend std::ostream& operator<<(std::ostream& os, const Data &p);

protected:

  boost::unordered_map< std::string, arma::mat> vector_data;
  
  boost::unordered_map< std::string, boost::spirit::hold_any> any_data;

  void make_copy(const Data &another);

};

//typedef Data Data;

//RCPP_EXPOSED_WRAP(Data);

#include <RcppArmadillo.h>

inline Data::Data()
{
  this->vector_data.clear();
  this->any_data.clear();
}

inline Data::~Data()
{

}

//Copy constructor for the Data class.
inline Data::Data(const Data &another)
{
  this->make_copy(another);
}

inline void Data::operator=(const Data &another)
{
  if(this == &another){ //if a==a
    return;
  }

  this->make_copy(another);
}

inline Data* Data::duplicate() const
{
  return( new Data(*this));
}

inline bool Data::is_empty() const
{
  if ( (this->vector_data.size()==0) && (this->any_data.size()==0) )
    return TRUE;
  else
    return FALSE;
}

/*
inline bool Data::operator==(const Data &another) const
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

inline bool Data::operator!=(const Data &another) const
{
  if (this->vector_data != another.vector_data)
  {
    return true;
  }
  
  if (this->any_data != another.any_data)
  {
    return true;
  }
  
  return false;
}
*/

inline arma::colvec Data::get_vector() const
{
  arma::colvec concatenated_vector;
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

inline arma::colvec Data::get_vector(const std::string &variable) const
{
  //arma::mat output = (*this)[variable];
  return arma::vectorise((*this)[variable]);
}

inline arma::colvec Data::get_vector(const std::vector<std::string> &variables) const
{
  arma::colvec concatenated_vector;
  for (auto i=variables.begin();
       i!=variables.end();
       ++i)
  {
    if (i==variables.begin())
      concatenated_vector = this->get_vector(*i);
    else
      concatenated_vector = join_rows(concatenated_vector,this->get_vector(*i));
  }
  return concatenated_vector;
}

inline void Data::make_copy(const Data &another)
{
  this->vector_data = another.vector_data;
  this->any_data = another.any_data;
}

inline arma::mat& Data::operator[](const std::string &variable)
{
  return(this->vector_data[variable]);
}

inline arma::mat Data::operator[](const std::string &variable) const
{
  boost::unordered_map< std::string, arma::mat>::const_iterator found =  this->vector_data.find(variable);

  if (found != this->vector_data.end())
    return(found->second);
  else
    return arma::mat();
}

inline arma::mat Data::operator[](const std::vector<std::string> &variables) const
{
  arma::mat concatenated_matrix;
  for (std::vector<std::string>::const_iterator i=variables.begin();
       i!=variables.end();
       ++i)
  {
    boost::unordered_map< std::string, arma::mat>::const_iterator found =  this->vector_data.find(*i);
    if (i==variables.begin())
    {
      if (found != this->vector_data.end())
        concatenated_matrix = found->second;
      else
        throw std::runtime_error("Data::operator[]: variable not found in Data.");
    }
    else
    {
      if (found != this->vector_data.end())
        concatenated_matrix = join_rows(concatenated_matrix,found->second);
      else
        throw std::runtime_error("Data::operator[]: variable not found in Data.");
    }
  }
  
  return concatenated_matrix;
}

inline boost::spirit::hold_any& Data::operator()(const std::string &variable)
{
  return(this->any_data[variable]);
}

inline boost::spirit::hold_any Data::operator()(const std::string &variable) const
{
  boost::unordered_map< std::string, boost::spirit::hold_any>::const_iterator found =  this->any_data.find(variable);

  if (found != this->any_data.end())
    return(found->second);
  else
    throw std::runtime_error("Data::operator(): variable not found in Data.");
}

inline Data Data::merge(const Data &another) const
{
  Data output = *this;
  
  output.vector_data.insert(another.vector_data.begin(), another.vector_data.end());
  
  output.any_data.insert(another.any_data.begin(), another.any_data.end());
  
  return output;
}

inline vector_parameter_iterator Data::vector_begin()
{
  return this->vector_data.begin();
}

inline vector_parameter_iterator Data::vector_end()
{
  return this->vector_data.end();
}

inline any_parameter_iterator Data::any_begin()
{
  return this->any_data.begin();
}

inline any_parameter_iterator Data::any_end()
{
  return this->any_data.end();
}

inline vector_parameter_const_iterator Data::vector_begin() const
{
  return this->vector_data.begin();
}

inline vector_parameter_const_iterator Data::vector_end() const
{
  return this->vector_data.end();
}

inline any_parameter_const_iterator Data::any_begin() const
{
  return this->any_data.begin();
}

inline any_parameter_const_iterator Data::any_end() const
{
  return this->any_data.end();
}

inline size_t Data::vector_size() const
{
  return this->vector_data.size();
}

inline size_t Data::any_size() const
{
  return this->any_data.size();
}

inline std::ostream& operator<<(std::ostream& os, const Data &p)
{
  boost::unordered_map< std::string, arma::mat>::const_iterator it;

  for (it=p.vector_data.begin();it!=p.vector_data.end();++it)
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

  for (it2=p.any_data.begin();it2!=p.any_data.end();++it2)
  {
    os << it2->first << ":";
    os << it2->second;
    os << std::endl;
  }

  return os;
}

#endif
