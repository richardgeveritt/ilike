// Include the class definition.
#include "mileometer.h"

// Include other class definitions that are needed.

// Public functions //

// Default constructor and destructor.

// Comment about function.
Mileometer::Mileometer(void)
{
	// Set all the pointers owned by this class to be something (possibly NULL).
}

// Comment about function.
Mileometer::~Mileometer(void)
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
Mileometer* Mileometer::duplicate(void) const
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
void Mileometer::operator=(const Mileometer &another)
{
	if(this==&another)//a==a
		return;
  
	this->current_index.clear();
	this->limits.clear();
	
  this->make_copy(another);
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
    if (*i==0)
    {
      return false;
    }
  }
  return true;
}

std::vector<double> Mileometer::get_current_values(const std::vector< std::vector<double> > &values) const
{
  if (this->current_index.size()!=values.size())
    throw std::runtime_error("Mileometer::get_current_value - index different length to values.");
  std::vector<double> current_value;
  current_value.reserve(values.size());
  for (size_t i=0; i<values.size(); ++i)
  {
    current_value.push_back(values[i][this->current_index[i]]);
  }
  return current_value;
}

// Increment current and return it.
void Mileometer::increment(void)
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
