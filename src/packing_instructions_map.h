//
//  packing_instructions_map.h
//  ilike_cpp
//
//  Created by Richard Everitt on 21/02/2023.
//

#ifndef packing_instructions_map_h
#define packing_instructions_map_h

#include "parameters.h"

namespace ilike
{
  /**
   * @file packing_instructions_map.h
   * @brief Defines the PackingInstructionsMap class.
   *
   * Provides packing instructions map functionality.
   *
   * @namespace ilike
   * @class PackingInstructionsMap
   * @brief The packing instructions map class.
   */


class PackingInstructionsMap
{
public:
  
  /**
   * @brief Default constructor for PackingInstructionsMap.
   */
  PackingInstructionsMap();
  /**
   * @brief Destructor for PackingInstructionsMap.
   */
  virtual ~PackingInstructionsMap();
  
  /**
   * @brief Copy constructor for PackingInstructionsMap.
   *
   * @param another The PackingInstructionsMap instance to copy from.
   */
  PackingInstructionsMap(const PackingInstructionsMap &another);
  /**
   * @brief Assignment operator for PackingInstructionsMap.
   *
   * @param another The PackingInstructionsMap instance to copy from.
   */
  PackingInstructionsMap& operator=(const PackingInstructionsMap &another);
  /**
   * @brief Copies the state of another PackingInstructionsMap into this object.
   *
   * @param another The PackingInstructionsMap instance to copy from.
   */
  void make_copy(const PackingInstructionsMap &another);
  
  void set_info(const Parameters &parameters,
                const std::vector<std::string> &variables);
  
  /**
   * @brief Performs the operator[] operation.
   *
   * @param variable The variable.
   *
   * @return The result.
   */
  std::pair<size_t,size_t>& operator[](const std::string &variable);
  /**
   * @brief Performs the operator[] operation.
   *
   * @param variable The variable.
   *
   * @return The result.
   */
  std::pair<size_t,size_t> operator[](const std::string &variable) const;
  
  boost::unordered_map< std::string, std::pair<size_t,size_t>> variable_indexed_start_and_end;
};
}

#endif /* packing_instructions_h */
