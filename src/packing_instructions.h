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
class PackingInstructions
{
public:
  
  PackingInstructions();
  virtual ~PackingInstructions();
  
  PackingInstructions(const PackingInstructions &another);
  PackingInstructions& operator=(const PackingInstructions &another);
  void make_copy(const PackingInstructions &another);
  
  void set_info(const Parameters &parameters,
                const std::vector<std::string> &variables);
  
  std::vector<std::string> states_names;
  std::vector<std::pair<size_t,size_t>> states_start_and_end;
};
}

#endif /* packing_instructions_h */
