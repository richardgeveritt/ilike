// Include the class definition.
#include "mileometer.h"

namespace ilike
{
// Include other class definitions that are needed.

// Public functions //

// Default constructor and destructor.

// Comment about function.
Mileometer::Mileometer()
{
  // Set all the pointers owned by this class to be something (possibly NULL).
}

// Comment about function.
Mileometer::~Mileometer()
{
  // Delete all the pointers owned by this class.
}

// Other constuctors.

// If the sizes vector is empty, then
Mileometer::Mileometer(const std::vector<size_t> &limitsin)
{
  this->limits = limitsin;
  this->current_index.reserve(this->limits.size());
  
  for (std::vector<size_t>::const_iterator j=this->limits.begin();
       j!=this->limits.end();
       ++j)
  {
    this->current_index.push_back(0);
  }
}

// Everything you need to copy the class.

// The copy constructor.
Mileometer::Mileometer(const Mileometer &another)
{
  this->make_copy(another);
}

// Returns a pointer to a base class of type Mileometer.
Mileometer* Mileometer::duplicate() const
{
  return new Mileometer(*this);
}

// Copy all the members of the class.
void Mileometer::make_copy(const Mileometer &another)
{
  // Copy all members, duplicating the memory where appropriate.
  this->current_index = another.current_index;
  this->limits = another.limits;
}

// The = operator.
Mileometer& Mileometer::operator=(const Mileometer &another)
{
  if(this==&another)//a==a
    return *this;
  
  this->current_index.clear();
  this->limits.clear();
  
  this->make_copy(another);
  
  return *this;
}

// The move constructor.
Mileometer::Mileometer(Mileometer &&another)
{
  this->make_copy(std::move(another));
}

// Copy all the members of the class.
void Mileometer::make_copy(Mileometer &&another)
{
  // Copy all members, duplicating the memory where appropriate.
  this->current_index = std::move(another.current_index);
  this->limits = std::move(another.limits);
  
  another.current_index = std::vector<size_t>();
  another.limits = std::vector<size_t>();
}

// The = operator.
Mileometer& Mileometer::operator=(Mileometer &&another)
{
  if(this==&another)//a==a
    return *this;
  
  this->current_index.clear();
  this->limits.clear();
  
  this->make_copy(std::move(another));
  
  return *this;
}

// Get the current index.
std::vector<size_t> Mileometer::get_current_index() const
{
  return this->current_index;
}

bool Mileometer::at_start() const
{
  for (std::vector<size_t>::const_iterator i=this->current_index.begin();
       i!=this->current_index.end();
       ++i)
  {
    if (*i!=0)
    {
      return false;
    }
  }
  return true;
}

std::vector<double> Mileometer::get_current_values(const std::vector< std::vector<double> > &values) const
{
  if (this->current_index.size()!=values.size())
    Rcpp::stop("Mileometer::get_current_value - index different length to values.");
  std::vector<double> current_value;
  current_value.reserve(values.size());
  for (size_t i=0; i<values.size(); ++i)
  {
    current_value.push_back(values[i][this->current_index[i]]);
  }
  return current_value;
}

size_t Mileometer::operator[](const size_t &i) const
{
  return this->current_index[i];
}

size_t Mileometer::size() const
{
  return this->current_index.size();
}

size_t Mileometer::back() const
{
  return this->current_index.back();
}

size_t Mileometer::get_previous_index() const
{
  return this->current_index.back();
}

size_t Mileometer::get_next_index() const
{
  if (this->current_index.back()>=this->limits.back()-1) // at end
    return 0;
  else
    return this->current_index.back()+1;
}

// Increment current and return it.
void Mileometer::increment()
{
  // Call increment on the end counter.
  this->increment(this->current_index.size()-1);
}

// Increment current dimension.
void Mileometer::increment(size_t dimension)
{
  // The current index.
  std::vector<size_t>::iterator the_current_index = current_index.begin()+dimension;
  
  // The current limit.
  std::vector<size_t>::iterator the_current_limit = limits.begin()+dimension;
  
  // If we pass the size of this dimension, also increment the previous dimension.
  // Unless we are already at the smallest dimension.
  if (*the_current_index != (*the_current_limit) - 1 )
  {
    // Increment.
    *the_current_index = *the_current_index + 1;
  }
  else
  {
    // Set back to the start.
    *the_current_index = 0;
    
    // Sort out the other dimensions.
    if (dimension!=0)
      this->increment(dimension-1);
  }
}

void Mileometer::reset_final_dimension(size_t new_number_of_states)
{
  this->limits.back() = new_number_of_states;
  this->current_index.back() = 0;
}

void Mileometer::reset()
{
  for (size_t i=0; i<this->current_index.size(); ++i)
  {
    current_index[i] = 0;
  }
}
}
