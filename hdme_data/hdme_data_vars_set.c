
#include "hdme_data.h"

void hdme_data_vars_set(char** vars, const char* name, slong k)
{
  strcpy(vars[k], name);
}
