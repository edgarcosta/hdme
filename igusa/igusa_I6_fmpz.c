
#include "igusa.h"

int igusa_I6_fmpz(fmpz_t I6, fmpz* I)
{
  fmpz_t I2, res, temp;
  int r;

  fmpz_init(I2);
  fmpz_init(res);
  fmpz_init(temp);

  r = igusa_I2_fmpz(I2, I);
  if (!r) return 0;
  
  /* Get I6 from I6prime */
  fmpz_mul_si(res, cov_I6prime(I), 2);
  fmpz_mul(temp, I2, cov_I4(I));
  fmpz_sub(res, res, temp);
  r = fmpz_divisible_si(res, 3);
  if (r)
    {
      fmpz_divexact_si(res, res, -3);
      fmpz_set(I6, res);
    }

  fmpz_clear(I2);
  fmpz_clear(res);
  fmpz_clear(temp);
  return r;
}
