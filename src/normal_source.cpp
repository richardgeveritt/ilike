#include "normal_source.h"

namespace ilike
{

thread_local NormalSource* installed_normal_source = NULL;
std::atomic<int> normal_sources_in_use(0);

}
