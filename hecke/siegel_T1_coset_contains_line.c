
#include "hecke.h"

int siegel_T1_coset_contains_line(const fmpz_mat_t m, const fmpz_mat_t L, slong ell)
{
  fmpz_mat_t r;
  nmod_mat_t red;
  int res;
  
  fmpz_mat_init(r, 4, 4);
  nmod_mat_init(red, 4, 4, ell);

  fmpz_mat_mul(r, m, L);
  fmpz_mat_get_nmod_mat(red, r);
  res = nmod_mat_is_zero(red);

  fmpz_mat_clear(r);
  nmod_mat_clear(red);
  return res;
}
