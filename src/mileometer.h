#ifndef MILEOMETER_H
#define MILEOMETER_H

// The point of this is to iterate through a vector of integers in the following way:
// 0 0 0 0
// 0 0 0 1
// 0 0 0 2
// 0 0 1 0
// 0 0 1 1
// 0 0 1 2
// 0 1 0 0
// 0 1 0 1
// etc.

#include <vector>
#include <stdexcept>

class Sequencer;
class EnsembleSequencer;

class Mileometer
{
public:

	Mileometer();
	virtual ~Mileometer();

	Mileometer(const std::vector<size_t> &limitsin);

	// Everything you need to copy the class.
	Mileometer(const Mileometer &another);
	Mileometer* duplicate() const;
	void make_copy(const Mileometer &another);
	void operator=(const Mileometer &another);

	// Get the current index.
	std::vector<size_t> get_current_index() const;
  
  std::vector<double> get_current_values(const std::vector< std::vector<double> > &values) const;
  
  bool at_start() const;

	// Increment current.
	void increment();

protected: // Things that can be accessed in this class and subclasses.

  friend Sequencer;
  friend EnsembleSequencer;
  
	std::vector<size_t> current_index;

	std::vector<size_t> limits;

	// Increment current dimension.
	void increment(size_t dimension);

private: // Things that can be accessed only by this class.

};

#endif