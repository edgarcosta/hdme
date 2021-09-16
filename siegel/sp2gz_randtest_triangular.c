
#include "siegel.h"

void
sp2gz_randtest_triangular(sp2gz_t m, flint_rand_t state, slong bits)
{
  fmpz_mat_t b, bt;
  bits = FLINT_MAX(bits, 1);
  
  fmpz_mat_init(b, m->g, m->g);
  fmpz_mat_init(bt, m->g, m->g);
  
  fmpz_mat_randbits(b, state, bits);
  fmpz_mat_transpose(bt, b);
  fmpz_mat_add(b, b, bt);
  fmpz_mat_scalar_tdiv_q_2exp(b, b, 1);
  sp2gz_one(m);
  fmpz_mat_set(&m->b, b);

  fmpz_mat_clear(b);
  fmpz_mat_clear(bt);
}
