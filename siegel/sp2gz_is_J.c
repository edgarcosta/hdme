
#include "siegel.h"

int sp2gz_is_J(const sp2gz_t m)
{
  fmpz_mat_t mc;
  int res;
  
  fmpz_mat_init(mc, m->g, m->g);
  fmpz_mat_neg(mc, &m->c);
  res = fmpz_mat_is_zero(&m->a)
    &&  fmpz_mat_is_one(&m->b)
    &&  fmpz_mat_is_one(mc)
    &&  fmpz_mat_is_zero(&m->d);

  fmpz_mat_clear(mc);
  return res;
}
