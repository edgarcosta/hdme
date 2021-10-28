
#include "siegel.h"

void fmpz_mat_J(fmpz_mat_t m)
{
  slong g = fmpz_mat_half_dim(m);
  fmpz_mat_t zero, one, minus_one;

  fmpz_mat_init(zero, g, g);
  fmpz_mat_init(one, g, g);
  fmpz_mat_init(minus_one, g, g);
  
  fmpz_mat_one(one);
  fmpz_mat_neg(minus_one, one);
  fmpz_mat_set_abcd(m, zero, one, minus_one, zero);
  
  fmpz_mat_clear(zero);
  fmpz_mat_clear(one);
  fmpz_mat_clear(minus_one);  
}
