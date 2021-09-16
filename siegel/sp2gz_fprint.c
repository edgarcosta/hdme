
#include "siegel.h"

void
sp2gz_fprint(FILE* file, const sp2gz_t m)
{
  fmpz_mat_t n;
  slong g = m->g;

  fmpz_mat_init(n, 2*g, 2*g);
  sp2gz_get_mat(n, m);
  fmpz_mat_fprint_pretty(file, n);
  fmpz_mat_clear(n);
}
