#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <RcppArmadillo.h>
//#include <RcppCommon.h>

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <boost/unordered_map.hpp>
#include <boost/any.hpp>
//#include <boost/spirit/home/support/detail/hold_any.hpp>

// was using boost::spirit::hold_any, but does not work with shared_ptr


// At some point, we might need a more flexible container for data and/or parameters.
// At this point, make another class that includes this one.

class Parameters;

typedef boost::unordered_map< std::string, std::pair<std::shared_ptr<arma::mat>,bool>>::iterator vector_parameter_iterator;
typedef boost::unordered_map< std::string, std::pair<std::shared_ptr<arma::mat>,bool>>::const_iterator vector_parameter_const_iterator;
typedef boost::unordered_map< std::string, std::pair<std::shared_ptr<boost::any>,bool>>::iterator any_parameter_iterator;
typedef boost::unordered_map< std::string, std::pair<std::shared_ptr<boost::any>,bool>>::const_iterator any_parameter_const_iterator;

RCPP_EXPOSED_CLASS(Parameters)

using boost::any_cast;

class Parameters
{
  
public:
  
  Parameters();
  Parameters(const std::string &variable_in,
             double value_in);
  Parameters(const std::string &variable_in,
             const arma::mat &value_in);
  //Parameters(const Rcpp::List &list_in);
  
  virtual ~Parameters();
  
  arma::mat& operator[](const std::string &variable);
  arma::mat operator[](const std::string &variable) const;
  arma::colvec operator[](const std::vector<std::string> &variables) const;
  
  //std::shared_ptr<arma::mat> get_vector_pointer(const std::string &variable);
  
  boost::any& operator()(const std::string &variable);
  boost::any operator()(const std::string &variable) const;
  
  //std::shared_ptr<boost::any> get_any_pointer(const std::string &variable);
  
  Parameters(const Parameters &another);
  Parameters& operator=(const Parameters &another);
  Parameters* duplicate() const;
  
  Parameters(Parameters &&another);
  Parameters& operator=(Parameters &&another);
  
  //bool operator==(const Parameters &another) const;
  //bool operator!=(const Parameters &another) const;
  
  bool is_empty() const;
  
  //arma::colvec get_vector() const;
  arma::colvec get_colvec(const std::string &variable) const;
  arma::colvec get_colvec(const std::vector<std::string> &variables) const;
  
  arma::rowvec get_rowvec(const std::string &variable) const;
  arma::rowvec get_rowvec(const std::vector<std::string> &variables) const;
  
  std::vector<arma::mat> get_matrices(const std::vector<std::string> &variables) const;
  
  Parameters row(size_t index) const;
  Parameters col(size_t index) const;
  
  Parameters rows(const arma::uvec &indices) const;
  Parameters cols(const arma::uvec &indices) const;
  
  size_t min_n_rows() const;
  size_t min_n_cols() const;
  
  Parameters merge(const Parameters &another) const;
  //void add_parameters(const Parameters &another);
  //void add_parameters_overwrite(const Parameters &another);
  
  void merge_with_fixed(const Parameters &conditioned_on_parameters);
  
  Parameters deep_copy() const;
  
  Parameters deep_copy_nonfixed() const;
  void self_deep_copy_nonfixed();
  
  //void overwrite_with_variables_in_argument(const Parameters &new_parameters);
  void deep_overwrite_with_variables_in_argument(const Parameters &new_parameters);
  
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
  
  std::vector<std::string> get_vector_variables() const;
  std::vector<std::string> get_any_variables() const;
  
  std::vector<std::string> get_nonfixed_vector_variables() const;
  std::vector<std::string> get_nonfixed_any_variables() const;
  
  std::vector<size_t> get_variable_n_elems(const std::vector<std::string> &variables) const;
  
  //Rcpp::List as_list() const;
  
  friend std::ostream& operator<<(std::ostream& os, const Parameters &p);
  
protected:
  
  boost::unordered_map< std::string, std::pair<std::shared_ptr<arma::mat>,bool>> vector_parameters;
  
  boost::unordered_map< std::string, std::pair<std::shared_ptr<boost::any>,bool>> any_parameters;
  
  void make_copy(const Parameters &another);
  void make_copy(Parameters &&another);
  
};

typedef Parameters Data;
typedef Parameters MatrixList;

//RCPP_EXPOSED_WRAP(Parameters);

#include <RcppArmadillo.h>

inline Parameters::Parameters()
{
  this->vector_parameters.clear();
  this->any_parameters.clear();
}

inline Parameters::Parameters(const std::string &variable_in,
                              double value_in)
{
  (*this)[variable_in] = value_in;
}

inline Parameters::Parameters(const std::string &variable_in,
                              const arma::mat &value_in)
{
  (*this)[variable_in] = value_in;
}

/*
inline Parameters::Parameters(const Rcpp::List &list_in)
{
  Rcpp::CharacterVector names = list_in.names();
  
  for (size_t i=0; i<names.size(); ++i)
  {
    for (size_t i=0; i<names.size(); ++i)
    {
      try
      {
        (*this)[Rcpp::as<std::string>(names[i])] = Rcpp::as<arma::mat>(list_in[i]);
      }
      catch (const std::runtime_error& re)
      {
        (*this)(Rcpp::as<std::string>(names[i])) = list_in[i];
      }
      catch(const std::exception& ex)
      {
        (*this)(Rcpp::as<std::string>(names[i]))= list_in[i];
      }
      catch(...)
      {
        std::cerr << "Parameters::Parameters(const Rcpp::List &list_in) - cannot read this element into Parameters." << std::endl;
      }
      
    }
  }
}
*/

inline Parameters::~Parameters()
{
  
}

//Copy constructor for the Parameters class.
inline Parameters::Parameters(const Parameters &another)
{
  this->make_copy(another);
}

inline Parameters& Parameters::operator=(const Parameters &another)
{
  if(this == &another){ //if a==a
    return *this;
  }
  
  this->make_copy(another);
  return *this;
}

//Move constructor for the Parameters class.
inline Parameters::Parameters(Parameters &&another)
{
  this->make_copy(std::move(another));
}

inline Parameters& Parameters::operator=(Parameters &&another)
{
  if(this == &another){ //if a==a
    return *this;
  }
  
  this->make_copy(std::move(another));
  return *this;
}

inline Parameters* Parameters::duplicate() const
{
  return( new Parameters(*this));
}

inline bool Parameters::is_empty() const
{
  if ( (this->vector_parameters.size()==0) && (this->any_parameters.size()==0) )
    return TRUE;
  else
    return FALSE;
}

/*
 inline bool Parameters::operator==(const Parameters &another) const
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
 
 inline bool Parameters::operator!=(const Parameters &another) const
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

/*
inline arma::colvec Parameters::get_vector() const
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
*/

inline arma::colvec Parameters::get_colvec(const std::string &variable) const
{
  //arma::mat output = (*this)[variable];
  return arma::vectorise((*this)[variable]);
}

inline arma::colvec Parameters::get_colvec(const std::vector<std::string> &variables) const
{
  arma::colvec concatenated_vector;
  for (auto i=variables.begin();
       i!=variables.end();
       ++i)
  {
    if (i==variables.begin())
      concatenated_vector = this->get_colvec(*i);
    else
      concatenated_vector = join_cols(concatenated_vector,this->get_colvec(*i));
  }
  return concatenated_vector;
}

inline arma::rowvec Parameters::get_rowvec(const std::string &variable) const
{
  //arma::mat output = (*this)[variable];
  return (*this)[variable].as_row();
}

inline arma::rowvec Parameters::get_rowvec(const std::vector<std::string> &variables) const
{
  arma::rowvec concatenated_vector;
  for (auto i=variables.begin();
       i!=variables.end();
       ++i)
  {
    if (i==variables.begin())
      concatenated_vector = this->get_rowvec(*i);
    else
      concatenated_vector = join_rows(concatenated_vector,this->get_rowvec(*i));
  }
  return concatenated_vector;
}

inline std::vector<arma::mat> Parameters::get_matrices(const std::vector<std::string> &variables) const
{
  std::vector<arma::mat> output;
  output.reserve(variables.size());
  for (auto i=variables.begin();
       i!=variables.end();
       ++i)
  {
    output.push_back((*this)[*i]);
  }
  return output;
}

inline Parameters Parameters::row(size_t index) const
{
  Parameters output;
  for (auto i=this->vector_begin(); i!=this->vector_end(); ++i)
  {
    output.vector_parameters.insert({i->first,std::pair<std::shared_ptr<arma::mat>,bool>(std::make_shared<arma::mat>(i->second.first->row(index)),i->second.second)});
  }
  return output;
}

inline Parameters Parameters::col(size_t index) const
{
  Parameters output;
  for (auto i=this->vector_begin(); i!=this->vector_end(); ++i)
  {
    output.vector_parameters.insert({i->first,std::pair<std::shared_ptr<arma::mat>,bool>(std::make_shared<arma::mat>(i->second.first->col(index)),i->second.second)});
  }
  return output;
}

inline Parameters Parameters::rows(const arma::uvec &indices) const
{
  Parameters output;
  for (auto i=this->vector_begin(); i!=this->vector_end(); ++i)
  {
    output.vector_parameters.insert({i->first,std::pair<std::shared_ptr<arma::mat>,bool>(std::make_shared<arma::mat>(i->second.first->rows(indices)),i->second.second)});
  }
  return output;
}

inline Parameters Parameters::cols(const arma::uvec &indices) const
{
  Parameters output;
  for (auto i=this->vector_begin(); i!=this->vector_end(); ++i)
  {
    output.vector_parameters.insert({i->first,std::pair<std::shared_ptr<arma::mat>,bool>(std::make_shared<arma::mat>(i->second.first->cols(indices)),i->second.second)});
  }
  return output;
}

inline size_t Parameters::min_n_rows() const
{
  size_t min = arma::datum::inf;
  for (auto i=this->vector_begin(); i!=this->vector_end(); ++i)
  {
    if (i->second.first->n_rows<min)
      min = i->second.first->n_rows;
  }
  return min;
}

inline size_t Parameters::min_n_cols() const
{
  size_t min = arma::datum::inf;
  for (auto i=this->vector_begin(); i!=this->vector_end(); ++i)
  {
    if (i->second.first->n_rows<min)
      min = i->second.first->n_cols;
  }
  return min;
}

inline void Parameters::make_copy(const Parameters &another)
{
  this->vector_parameters = another.vector_parameters;
  this->any_parameters = another.any_parameters;
}

inline void Parameters::make_copy(Parameters &&another)
{
  this->vector_parameters = std::move(another.vector_parameters);
  this->any_parameters = std::move(another.any_parameters);
  
  another.vector_parameters = boost::unordered_map< std::string, std::pair<std::shared_ptr<arma::mat>,bool>>();
  another.any_parameters = boost::unordered_map< std::string, std::pair<std::shared_ptr<boost::any>,bool>>();
}

inline arma::mat& Parameters::operator[](const std::string &variable)
{
  
  // First try.
  
  /*
  auto found = this->vector_parameters.find(variable);
  
  // if doesn't yet exist, make shared pointer, and set to be not fixed
  if (found==this->vector_parameters.end())
  {
    std::pair<std::shared_ptr<arma::mat>,bool>& new_element = this->vector_parameters["variable"];
    new_element = std::pair<std::shared_ptr<arma::mat>,bool>(std::make_shared<arma::mat>(),false);
    return *new_element.first;
  }
  else
  {
    if (found->second.second==false)
    {
      // if does exist, and is not fixed, just pass reference to memory that the shard pointer points to
      return *found->second.first;
    }
    else
    {
      // if does exist, and is fixed, throw error
      Rcpp::stop("Parameters::operator[] - parameter is fixed and cannot be changed.");
    }
  }
  */
  
  // Second try. 18
  
  
  auto found = this->vector_parameters.find(variable);
  
  // if doesn't yet exist, make shared pointer, and set to be not fixed
  if (found==this->vector_parameters.end())
  {
    return *this->vector_parameters.insert({variable,std::pair<std::shared_ptr<arma::mat>,bool>(std::make_shared<arma::mat>(),false)}).first->second.first;
    
    //std::pair<std::shared_ptr<arma::mat>,bool>& new_element = this->vector_parameters["variable"];
    //new_element = std::pair<std::shared_ptr<arma::mat>,bool>(std::make_shared<arma::mat>(),false);
    //return *new_element.first;
  }
  else
  {
    if (found->second.second==false)
    {
      // if does exist, and is not fixed, just pass reference to memory that the shard pointer points to
      return *found->second.first;
    }
    else
    {
      // if does exist, and is fixed, throw error
      Rcpp::stop("Parameters::operator[] - parameter is fixed and cannot be changed.");
    }
  }
  
  
  // Third try. 22.1
  
  /*
  auto new_element_result = this->vector_parameters.insert({variable,std::pair<std::shared_ptr<arma::mat>,bool>(std::make_shared<arma::mat>(),false)});
  
  // if doesn't yet exist, make shared pointer, and set to be not fixed
  if (new_element_result.second)
  {
    return *new_element_result.first->second.first;// .first->second.first;
    
    //std::pair<std::shared_ptr<arma::mat>,bool>& new_element = this->vector_parameters["variable"];
    //new_element = std::pair<std::shared_ptr<arma::mat>,bool>(std::make_shared<arma::mat>(),false);
    //return *new_element.first;
  }
  else
  {
    auto found = this->vector_parameters.find(variable);
    if (found->second.second==false)
    {
      // if does exist, and is not fixed, just pass reference to memory that the shard pointer points to
      return *found->second.first;
    }
    else
    {
      // if does exist, and is fixed, throw error
      Rcpp::stop("Parameters::operator[] - parameter is fixed and cannot be changed.");
    }
  }
  */
  
  // Fourth try. 19.0
  /*
  auto found = this->vector_parameters.find(variable);

  // if doesn't yet exist, make shared pointer, and set to be not fixed
  if (found==this->vector_parameters.end())
  {
    return *this->vector_parameters.emplace(variable,std::pair<std::shared_ptr<arma::mat>,bool>(std::make_shared<arma::mat>(),false)).first->second.first;

  //std::pair<std::shared_ptr<arma::mat>,bool>& new_element = this->vector_parameters["variable"];
  //new_element = std::pair<std::shared_ptr<arma::mat>,bool>(std::make_shared<arma::mat>(),false);
  //return *new_element.first;
  }
  else
  {
    if (found->second.second==false)
    {
      // if does exist, and is not fixed, just pass reference to memory that the shard pointer points to
      return *found->second.first;
    }
    else
    {
      // if does exist, and is fixed, throw error
      Rcpp::stop("Parameters::operator[] - parameter is fixed and cannot be changed.");
    }
  }
  */
  
  // Fifth try. 22.6
  
  /*
  auto new_element_result = this->vector_parameters.emplace(variable,std::pair<std::shared_ptr<arma::mat>,bool>(std::make_shared<arma::mat>(),false));
  
  // if doesn't yet exist, make shared pointer, and set to be not fixed
  if (new_element_result.second)
  {
    return *new_element_result.first->second.first;// .first->second.first;
    
    //std::pair<std::shared_ptr<arma::mat>,bool>& new_element = this->vector_parameters["variable"];
    //new_element = std::pair<std::shared_ptr<arma::mat>,bool>(std::make_shared<arma::mat>(),false);
    //return *new_element.first;
  }
  else
  {
    auto found = this->vector_parameters.find(variable);
    if (found->second.second==false)
    {
      // if does exist, and is not fixed, just pass reference to memory that the shard pointer points to
      return *found->second.first;
    }
    else
    {
      // if does exist, and is fixed, throw error
      Rcpp::stop("Parameters::operator[] - parameter is fixed and cannot be changed.");
    }
  }
  */
    
}

inline arma::mat Parameters::operator[](const std::string &variable) const
{
  auto found = this->vector_parameters.find(variable);
  if (found != this->vector_parameters.end())
    return(*found->second.first);
  else
    Rcpp::stop("Parameters::operator[]: variable '" + variable + "' not found in Parameters.");
}

inline arma::colvec Parameters::operator[](const std::vector<std::string> &variables) const
{
  arma::mat concatenated_matrix;
  for (std::vector<std::string>::const_iterator i=variables.begin();
       i!=variables.end();
       ++i)
  {
    auto found = this->vector_parameters.find(*i);
    if (i==variables.begin())
    {
      if (found != this->vector_parameters.end())
        concatenated_matrix = arma::vectorise(*found->second.first);
      else
        Rcpp::stop("Parameters::operator[]: variable not found in Parameters.");
    }
    else
    {
      if (found != this->vector_parameters.end())
        concatenated_matrix = join_cols(concatenated_matrix,arma::vectorise(*found->second.first));
      else
        Rcpp::stop("Parameters::operator[]: variable not found in Parameters.");
    }
  }
  
  return concatenated_matrix;
}

/*
inline std::shared_ptr<arma::mat> Parameters::get_vector_pointer(const std::string &variable)
{
  auto found = this->vector_parameters.find(variable);
  
  // if doesn't yet exist, make shared pointer, and set to be not fixed
  if (found==this->vector_parameters.end())
  {
    return this->vector_parameters.insert({variable,std::pair<std::shared_ptr<arma::mat>,bool>(std::make_shared<arma::mat>(),false)}).first->second.first;
    
    //std::pair<std::shared_ptr<arma::mat>,bool>& new_element = this->vector_parameters["variable"];
    //new_element = std::pair<std::shared_ptr<arma::mat>,bool>(std::make_shared<arma::mat>(),false);
    //return *new_element.first;
  }
  else
  {
    if (found->second.second==false)
    {
      // if does exist, and is not fixed, just pass reference to memory that the shard pointer points to
      return found->second.first;
    }
    else
    {
      // if does exist, and is fixed, throw error
      Rcpp::stop("Parameters::operator[] - parameter is fixed and cannot be changed.");
    }
  }
}
*/

inline boost::any& Parameters::operator()(const std::string &variable)
{
  auto found = this->any_parameters.find(variable);
  
  // if doesn't yet exist, make shared pointer, and set to be not fixed
  if (found==this->any_parameters.end())
  {
    return *this->any_parameters.insert({variable,std::pair<std::shared_ptr<boost::any>,bool>(std::make_shared<boost::any>(),false)}).first->second.first;
    
    //std::pair<std::shared_ptr<arma::mat>,bool>& new_element = this->vector_parameters["variable"];
    //new_element = std::pair<std::shared_ptr<arma::mat>,bool>(std::make_shared<arma::mat>(),false);
    //return *new_element.first;
  }
  else
  {
    if (found->second.second==false)
    {
      // if does exist, and is not fixed, just pass reference to memory that the shard pointer points to
      return *found->second.first;
    }
    else
    {
      // if does exist, and is fixed, throw error
      Rcpp::stop("Parameters::operator() - parameter is fixed and cannot be changed.");
    }
  }
}

inline boost::any Parameters::operator()(const std::string &variable) const
{
  auto found = this->any_parameters.find(variable);
  
  if (found != this->any_parameters.end())
    return(*found->second.first);
  else
    Rcpp::stop("Parameters::operator(): variable not found in Parameters.");
}

/*
inline std::shared_ptr<boost::any> Parameters::get_any_pointer(const std::string &variable)
{
  auto new_element_result = this->any_parameters.emplace(variable,std::pair<std::shared_ptr<boost::any>,bool>(std::make_shared<boost::any>(),false));
  
  // if doesn't yet exist, make shared pointer, and set to be not fixed
  if (new_element_result.second)
  {
    return new_element_result.first->second.first;// .first->second.first;
    
    //std::pair<std::shared_ptr<arma::mat>,bool>& new_element = this->vector_parameters["variable"];
    //new_element = std::pair<std::shared_ptr<arma::mat>,bool>(std::make_shared<arma::mat>(),false);
    //return *new_element.first;
  }
  else
  {
    auto found = this->any_parameters.find(variable);
    if (found->second.second==false)
    {
      // if does exist, and is not fixed, just pass reference to memory that the shard pointer points to
      return found->second.first;
    }
    else
    {
      // if does exist, and is fixed, throw error
      Rcpp::stop("Parameters::operator[] - parameter is fixed and cannot be changed.");
      Rcpp::stop("Parameters::operator[] - parameter is fixed and cannot be changed.");
    }
  }
}
*/

inline Parameters Parameters::merge(const Parameters &another) const
{
  Parameters output;
  output.vector_parameters = this->vector_parameters;
  output.any_parameters = this->any_parameters;
  
  output.vector_parameters.insert(another.vector_parameters.begin(), another.vector_parameters.end());
  
  output.any_parameters.insert(another.any_parameters.begin(), another.any_parameters.end());
  
  return output;
}

inline void Parameters::merge_with_fixed(const Parameters &conditioned_on_parameters)
{
  for (auto i=conditioned_on_parameters.vector_begin(); i!=conditioned_on_parameters.vector_end(); ++i)
  {
    this->vector_parameters[i->first] = std::pair<std::shared_ptr<arma::mat>,bool>(i->second.first,true);
  }
  
  for (auto i=conditioned_on_parameters.any_begin(); i!=conditioned_on_parameters.any_end(); ++i)
  {
    this->any_parameters[i->first] = std::pair<std::shared_ptr<boost::any>,bool>(i->second.first,true);
  }
}

inline Parameters Parameters::deep_copy() const
{
  Parameters output;
  
  for (auto i=this->vector_begin();
       i!=this->vector_end();
       ++i)
  {
    if (i->second.second==true) // fixed
    {
      output.vector_parameters.insert({i->first,std::pair<std::shared_ptr<arma::mat>,bool>(std::make_shared<arma::mat>(*i->second.first),true)});
    }
    else // non-fixed
    {
      output.vector_parameters.insert({i->first,std::pair<std::shared_ptr<arma::mat>,bool>(std::make_shared<arma::mat>(*i->second.first),false)});
    }
  }
  
  for (auto i=this->any_begin();
       i!=this->any_end();
       ++i)
  {
    if (i->second.second==true) // fixed
    {
      output.any_parameters.insert({i->first,std::pair<std::shared_ptr<boost::any>,bool>(std::make_shared<boost::any>(*i->second.first),true)});
    }
    else // non-fixed
    {
      output.any_parameters.insert({i->first,std::pair<std::shared_ptr<boost::any>,bool>(std::make_shared<boost::any>(*i->second.first),false)});
    }
  }
  
  return output;
}

inline Parameters Parameters::deep_copy_nonfixed() const
{
  Parameters output;
  
  for (auto i=this->vector_begin();
       i!=this->vector_end();
       ++i)
  {
    if (i->second.second==true) // fixed, shallow copy
    {
      output.vector_parameters.insert({i->first,std::pair<std::shared_ptr<arma::mat>,bool>(i->second.first,true)});
    }
    else // non-fixed, deep copy
    {
      output.vector_parameters.insert({i->first,std::pair<std::shared_ptr<arma::mat>,bool>(std::make_shared<arma::mat>(*i->second.first),false)});
    }
  }
  
  for (auto i=this->any_begin();
       i!=this->any_end();
       ++i)
  {
    if (i->second.second==true) // fixed, shallow copy
    {
      output.any_parameters.insert({i->first,std::pair<std::shared_ptr<boost::any>,bool>(i->second.first,true)});
    }
    else // non-fixed, deep copy
    {
      output.any_parameters.insert({i->first,std::pair<std::shared_ptr<boost::any>,bool>(std::make_shared<boost::any>(*i->second.first),false)});
    }
  }
  
  return output;
}

inline void Parameters::self_deep_copy_nonfixed()
{
  
  for (auto i=this->vector_begin();
       i!=this->vector_end();
       ++i)
  {
    if (i->second.second==false) // not fixed, deep copy
    {
      i->second.first = std::make_shared<arma::mat>(*i->second.first);
    }
  }
  
  for (auto i=this->any_begin();
       i!=this->any_end();
       ++i)
  {
    if (i->second.second==false) // not fixed, deep copy
    {
      i->second.first = std::make_shared<boost::any>(*i->second.first);
    }
  }
  
}

/*
inline void Parameters::overwrite_with_variables_in_argument(const Parameters &new_parameters)
{
  for (auto i=new_parameters.vector_begin(); i!=new_parameters.vector_end(); ++i)
  {
    this->get_vector_pointer(i->first) = i->second.first; // need similar function that takes ptr, not mat
  }
  
  for (auto i=new_parameters.any_begin(); i!=new_parameters.any_end(); ++i)
  {
    this->get_any_pointer(i->first) = i->second.first;
  }
}
*/

inline void Parameters::deep_overwrite_with_variables_in_argument(const Parameters &new_parameters)
{
  for (auto i=new_parameters.vector_begin(); i!=new_parameters.vector_end(); ++i)
  {
    auto found = this->vector_parameters.find(i->first);
    
    // if doesn't yet exist, make shared pointer, and set to be not fixed
    if (found==this->vector_parameters.end())
    {
      this->vector_parameters.insert({i->first,std::pair<std::shared_ptr<arma::mat>,bool>(std::make_shared<arma::mat>(*i->second.first),false)});
      
      //std::pair<std::shared_ptr<arma::mat>,bool>& new_element = this->vector_parameters["variable"];
      //new_element = std::pair<std::shared_ptr<arma::mat>,bool>(std::make_shared<arma::mat>(),false);
      //return *new_element.first;
    }
    else
    {
      if (found->second.second==false)
      {
        // if does exist, and is not fixed, just pass reference to memory that the shard pointer points to
        found->second.first = std::make_shared<arma::mat>(*i->second.first);
      }
      else
      {
        // if does exist, and is fixed, throw error
        Rcpp::stop("Parameters::deep_overwrite_with_variables_in_argument - parameter is fixed and cannot be changed.");
      }
    }
    
    //this->get_vector_pointer(i->first) = std::make_shared<arma::mat>(*i->second.first);
    //(*this)[i->first] = std::make_shared<arma::mat>(*i->second.first);
  }
  
  for (auto i=new_parameters.any_begin(); i!=new_parameters.any_end(); ++i)
  {
    auto found = this->any_parameters.find(i->first);
    
    // if doesn't yet exist, make shared pointer, and set to be not fixed
    if (found==this->any_parameters.end())
    {
      this->any_parameters.insert({i->first,std::pair<std::shared_ptr<boost::any>,bool>(std::make_shared<boost::any>(*i->second.first),false)});
      
      //std::pair<std::shared_ptr<arma::mat>,bool>& new_element = this->vector_parameters["variable"];
      //new_element = std::pair<std::shared_ptr<arma::mat>,bool>(std::make_shared<arma::mat>(),false);
      //return *new_element.first;
    }
    else
    {
      if (found->second.second==false)
      {
        // if does exist, and is not fixed, just pass reference to memory that the shard pointer points to
        found->second.first = std::make_shared<boost::any>(*i->second.first);
      }
      else
      {
        // if does exist, and is fixed, throw error
        Rcpp::stop("Parameters::deep_overwrite_with_variables_in_argument - parameter is fixed and cannot be changed.");
      }
    }
    
    //this->get_vector_pointer(i->first) = std::make_shared<arma::mat>(*i->second.first);
    //(*this)[i->first] = std::make_shared<arma::mat>(*i->second.first);
  }
}

/*
inline void Parameters::add_parameters(const Parameters &another)
{
  this->vector_parameters.insert(another.vector_parameters.begin(), another.vector_parameters.end());
  this->any_parameters.insert(another.any_parameters.begin(), another.any_parameters.end());
}

inline void Parameters::add_parameters_overwrite(const Parameters &another)
{
  for (auto i=another.vector_begin(); i!=another.vector_end(); ++i)
  {
    this->vector_parameters[i->first] = i->second;
  }
  
  for (auto i=another.any_begin(); i!=another.any_end(); ++i)
  {
    this->any_parameters[i->first] = i->second;
  }
}
*/

inline vector_parameter_iterator Parameters::vector_begin()
{
  return this->vector_parameters.begin();
}

inline vector_parameter_iterator Parameters::vector_end()
{
  return this->vector_parameters.end();
}

inline any_parameter_iterator Parameters::any_begin()
{
  return this->any_parameters.begin();
}

inline any_parameter_iterator Parameters::any_end()
{
  return this->any_parameters.end();
}

inline vector_parameter_const_iterator Parameters::vector_begin() const
{
  return this->vector_parameters.begin();
}

inline vector_parameter_const_iterator Parameters::vector_end() const
{
  return this->vector_parameters.end();
}

inline any_parameter_const_iterator Parameters::any_begin() const
{
  return this->any_parameters.begin();
}

inline any_parameter_const_iterator Parameters::any_end() const
{
  return this->any_parameters.end();
}

inline size_t Parameters::vector_size() const
{
  return this->vector_parameters.size();
}

inline size_t Parameters::any_size() const
{
  return this->any_parameters.size();
}

inline std::vector<std::string> Parameters::get_vector_variables() const
{
  std::vector<std::string> variables;
  for (auto it=this->vector_parameters.begin();it!=this->vector_parameters.end();++it)
  {
    variables.push_back(it->first);
  }
  return variables;
}

inline std::vector<std::string> Parameters::get_any_variables() const
{
  std::vector<std::string> variables;
  for (auto it=this->any_parameters.begin();it!=this->any_parameters.end();++it)
  {
    variables.push_back(it->first);
  }
  return variables;
}

inline std::vector<std::string> Parameters::get_nonfixed_vector_variables() const
{
  std::vector<std::string> variables;
  for (auto it=this->vector_parameters.begin();it!=this->vector_parameters.end();++it)
  {
    if (it->second.second==false)
      variables.push_back(it->first);
  }
  return variables;
}

inline std::vector<std::string> Parameters::get_nonfixed_any_variables() const
{
  std::vector<std::string> variables;
  for (auto it=this->any_parameters.begin();it!=this->any_parameters.end();++it)
  {
    if (it->second.second==false)
      variables.push_back(it->first);
  }
  return variables;
}

inline std::vector<size_t> Parameters::get_variable_n_elems(const std::vector<std::string> &variables) const
{
  std::vector<size_t> n_elems;
  n_elems.reserve(variables.size());
  for (size_t i=0; i<variables.size(); ++i)
  {
    n_elems.push_back((*this)[variables[i]].n_elem);
  }
  return n_elems;
}

/*
inline Rcpp::List Parameters::as_list() const
{
  Rcpp::List output;
  
  for (auto it=this->vector_parameters.begin();it!=this->vector_parameters.end();++it)
  {
    output[it->first] = *it->second.first;
  }
  
  
  for (auto it2=this->any_parameters.begin();it2!=this->any_parameters.end();++it2)
  {
    output[it2->first] = boost::any_cast<SEXP>(*it2->second.first);
  }
  
  
  return output;
}
*/

inline Parameters pow(const Parameters &p, double power)
{
  Parameters output;
  for (auto it=p.vector_begin();it!=p.vector_end();++it)
  {
    output[it->first] = arma::pow(p[it->first],power);
  }
  
  return output;
}

inline Parameters operator*(double scale, const Parameters &p)
{
  Parameters output = p;
  for (auto it=p.vector_begin();it!=p.vector_end();++it)
  {
    output[it->first] = scale*p[it->first];
  }
  return output;
}

inline std::ostream& operator<<(std::ostream& os, const Parameters &p)
{
  for (auto it=p.vector_parameters.begin();it!=p.vector_parameters.end();++it)
  {
    if (it->first!="")
    {
      os << it->first << "=";
      for (arma::mat::const_iterator i=it->second.first->begin(); i!=it->second.first->end(); ++i)
      {
        if (i==it->second.first->begin())
          os << *i;
        else
          os << ", " << *i;
      }
      os << ";";
    }
  }

  for (auto it2=p.any_parameters.begin();it2!=p.any_parameters.end();++it2)
  {
    if (it2->first!="")
    {
      os << it2->first << "=";
      os << it2->second.first;
      os << std::endl;
    }
  }
  
  return os;
}

#endif
