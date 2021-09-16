
#include "igusa.h"

void igusa_from_cov_fmpz(fmpq* j, const fmpz* I)
{
  fmpz_t num, den;
  fmpz_init(num);
  fmpz_init(den);

  fmpz_mul(num, &I[1], &I[2]);
  fmpz_set(den, &I[3]);
  fmpq_set_fmpz_frac(&j[0], num, den);

  fmpz_mul(num, &I[0], &I[1]);
  fmpz_mul(num, num, &I[1]);
  fmpq_set_fmpz_frac(&j[1], num, den);

  fmpz_pow_ui(num, &I[1], 5);
  fmpz_pow_ui(den, &I[3], 2);
  fmpq_set_fmpz_frac(&j[2], num, den);

  fmpz_clear(num);
  fmpz_clear(den);
}
