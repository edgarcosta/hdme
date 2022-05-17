
#include "hecke.h"

slong siegel_nb_T1_cosets(slong p)
{
  return p + n_pow(p, 2) + n_pow(p, 3) + n_pow(p, 4);
}
