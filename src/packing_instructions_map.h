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
class PackingInstructionsMap
{
public:
  
  PackingInstructionsMap();
  virtual ~PackingInstructionsMap();
  
  PackingInstructionsMap(const PackingInstructionsMap &another);
  PackingInstructionsMap& operator=(const PackingInstructionsMap &another);
  void make_copy(const PackingInstructionsMap &another);
  
  void set_info(const Parameters &parameters,
                const std::vector<std::string> &variables);
  
  std::pair<size_t,size_t>& operator[](const std::string &variable);
  std::pair<size_t,size_t> operator[](const std::string &variable) const;
  
  boost::unordered_map< std::string, std::pair<size_t,size_t>> variable_indexed_start_and_end;
};
}

#endif /* packing_instructions_h */
