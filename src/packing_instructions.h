//
//  packing_instructions.h
//  ilike_cpp
//
//  Created by Richard Everitt on 21/02/2023.
//

#ifndef packing_instructions_h
#define packing_instructions_h

#include "parameters.h"

namespace ilike
{
  /**
   * @file packing_instructions.h
   * @brief Defines the PackingInstructions class.
   *
   * Provides packing instructions functionality.
   *
   * @namespace ilike
   * @class PackingInstructions
   * @brief The packing instructions class.
   */


class PackingInstructions
{
public:
  
  /**
   * @brief Default constructor for PackingInstructions.
   */
  PackingInstructions();
  /**
   * @brief Destructor for PackingInstructions.
   */
  virtual ~PackingInstructions();
  
  /**
   * @brief Copy constructor for PackingInstructions.
   *
   * @param another The PackingInstructions instance to copy from.
   */
  PackingInstructions(const PackingInstructions &another);
  /**
   * @brief Assignment operator for PackingInstructions.
   *
   * @param another The PackingInstructions instance to copy from.
   */
  PackingInstructions& operator=(const PackingInstructions &another);
  /**
   * @brief Copies the state of another PackingInstructions into this object.
   *
   * @param another The PackingInstructions instance to copy from.
   */
  void make_copy(const PackingInstructions &another);
  
  void set_info(const Parameters &parameters,
                const std::vector<std::string> &variables);
  
  std::vector<std::string> states_names;
  std::vector<std::pair<size_t,size_t>> states_start_and_end;
};
}

#endif /* packing_instructions_h */
