#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <RcppArmadilloForward.h>
#include <RcppCommon.h>

#include <iostream>
#include <vector>
#include <boost/unordered_map.hpp>


// At some point, we might need a more flexible container for data and/or parameters.
// At this point, make another class that includes this one.

class Parameters;

RCPP_EXPOSED_CLASS(Parameters);

class Parameters
{

public:

  Parameters();

  virtual ~Parameters();

  arma::colvec& operator[](const std::string &variable);
  arma::colvec operator[](const std::string &variable) const;

  Parameters(const Parameters &another);
  void operator=(const Parameters &another);
  Parameters* duplicate() const;

  friend std::ostream& operator<<(std::ostream& os, const Parameters &p);

protected:

  boost::unordered_map< std::string, arma::colvec> parameters;

  void make_copy(const Parameters &another);

};

typedef Parameters Data;

//RCPP_EXPOSED_WRAP(Parameters);

#include <RcppArmadillo.h>

inline Parameters::Parameters()
{
  this->parameters.clear();
}

inline Parameters::~Parameters()
{

}

//Copy constructor for the Parameters class.
inline Parameters::Parameters(const Parameters &another)
{
  this->make_copy(another);
}

inline void Parameters::operator=(const Parameters &another)
{
  if(this == &another){ //if a==a
    return;
  }

  this->make_copy(another);
}

inline Parameters* Parameters::duplicate(void)const
{
  return( new Parameters(*this));
}

inline void Parameters::make_copy(const Parameters &another)
{
  this->parameters = another.parameters;
}

inline arma::colvec& Parameters::operator[](const std::string &variable)
{
  return(this->parameters[variable]);
}

inline arma::colvec Parameters::operator[](const std::string &variable) const
{
  boost::unordered_map< std::string, arma::colvec>::const_iterator found =  this->parameters.find(variable);

  if (found != this->parameters.end())
    return(found->second);
  else
    throw std::runtime_error("Parameters::operator[]: variable not found in Parameters.");
}

inline std::ostream& operator<<(std::ostream& os, const Parameters &p)
{
  boost::unordered_map< std::string, arma::colvec>::const_iterator it;

  for (it=p.parameters.begin();it!=p.parameters.end();++it)
  {
    os << it->first << ":";
    for (arma::colvec::const_iterator i=it->second.begin(); i!=it->second.end(); ++i)
    {
      if (i==it->second.begin())
        os << *i;
      else
        os << "," << *i;
    }
    os << " ";
  }

  return os;
}

#endif
