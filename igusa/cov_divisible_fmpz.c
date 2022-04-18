
#include "igusa.h"

int cov_divisible_fmpz(fmpz* I, const fmpz_t scal)
{
  fmpz_t f;
  int r;

  fmpz_init(f);

  r = fmpz_divisible(&I[0], scal);
  fmpz_pow_ui(f, scal, 2);
  r = r && fmpz_divisible(&I[1], f);
  fmpz_pow_ui(f, scal, 3);
  r = r && fmpz_divisible(&I[2], f);
  fmpz_pow_ui(f, scal, 5);
  r = r && fmpz_divisible(&I[3], f);

  fmpz_clear(f);
  return r;
}
