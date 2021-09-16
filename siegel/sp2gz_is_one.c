
#include "siegel.h"

int
sp2gz_is_one(const sp2gz_t m)
{
  return fmpz_mat_is_one(&m->a)
    &&   fmpz_mat_is_zero(&m->b)
    &&   fmpz_mat_is_zero(&m->c)
    &&   fmpz_mat_is_one(&m->d);
}
