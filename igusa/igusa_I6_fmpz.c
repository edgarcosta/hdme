
#include "igusa.h"

int igusa_I6_fmpz(fmpz_t I6, fmpz* I)
{
  fmpz_t res, temp;
  int r;

  fmpz_init(res);
  fmpz_init(temp);
  
  /* Get I6 from I6prime */
  fmpz_mul_si(res, &I[2], 2);
  fmpz_mul(temp, &I[0], &I[1]);
  fmpz_sub(res, res, temp);
  r = fmpz_divisible_si(res, 3);
  if (r)
    {
      fmpz_divexact_si(res, res, -3);
      fmpz_set(I6, res);
    }

  fmpz_clear(res);
  fmpz_clear(temp);
  return r;
}
