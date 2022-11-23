
#include "siegel.h"

int fmpz_mat_is_symplectic(const fmpz_mat_t m)
{
  slong g = fmpz_mat_half_dim(m);
  fmpz_mat_t a, b, c, d;
  fmpz_mat_t prod1, prod2;
  int res;

  fmpz_mat_init(a, g, g);
  fmpz_mat_init(b, g, g);
  fmpz_mat_init(c, g, g);
  fmpz_mat_init(d, g, g);
  fmpz_mat_init(prod1, g, g);
  fmpz_mat_init(prod2, g, g);

  fmpz_mat_get_a(a, m);
  fmpz_mat_get_b(b, m);
  fmpz_mat_get_c(c, m);
  fmpz_mat_get_d(d, m);

  fmpz_mat_transpose(prod1, a);
  fmpz_mat_mul(prod1, prod1, c);
  fmpz_mat_transpose(prod2, c);
  fmpz_mat_mul(prod2, prod2, a);
  fmpz_mat_sub(prod1, prod1, prod2);
  res = fmpz_mat_is_zero(prod1);

  fmpz_mat_transpose(prod1, b);
  fmpz_mat_mul(prod1, prod1, d);
  fmpz_mat_transpose(prod2, d);
  fmpz_mat_mul(prod2, prod2, b);
  fmpz_mat_sub(prod1, prod1, prod2);
  res = res && fmpz_mat_is_zero(prod1);

  fmpz_mat_transpose(prod1, a);
  fmpz_mat_mul(prod1, prod1, d);
  fmpz_mat_transpose(prod2, c);
  fmpz_mat_mul(prod2, prod2, b);
  fmpz_mat_sub(prod1, prod1, prod2);
  res = res && fmpz_mat_is_one(prod1);

  fmpz_mat_clear(a);
  fmpz_mat_clear(b);
  fmpz_mat_clear(c);
  fmpz_mat_clear(d);
  fmpz_mat_clear(prod1);
  fmpz_mat_clear(prod2);

  return res;
}
