
#include "siegel.h"

void sp2gz_set_diagonal(sp2gz_t m, const fmpz_mat_t u)
{
  fmpz_t den;
  fmpz_init(den);

  sp2gz_one(m);
  fmpz_mat_set(&m->a, u);
  fmpz_mat_inv(&m->d, den, &m->a);
  fmpz_mat_transpose(&m->d, &m->d);
  
  if (!fmpz_is_one(den))
    {
      fmpz_neg(den, den);
      fmpz_mat_neg(&m->d, &m->d);
    }
  if (!fmpz_is_one(den))
    {
      flint_fprintf(stderr, "Non-invertible matrix in sp2gz_set_diagonal:\n");
      fmpz_mat_fprint_pretty(stderr, u);
      flint_fprintf(stderr, "\n");
      flint_abort();
    }
  fmpz_clear(den);
}
