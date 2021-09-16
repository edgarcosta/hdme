
#include "siegel.h"

void
sp2gz_randtest_diagonal(sp2gz_t m, flint_rand_t state, slong bits)
{
  fmpz_t den;
  fmpz_init(den);
  
  bits = FLINT_MAX(bits, 1);
  sp2gz_one(m);
  fmpz_mat_randops(&m->a, state, 2 * bits * m->g);
  fmpz_mat_inv(&m->d, den, &m->a);
  fmpz_mat_transpose(&m->d, &m->d);
  if (!fmpz_is_one(den))
    fmpz_mat_neg(&m->d, &m->d);

  fmpz_clear(den);
}
