
#include "igusa.h"

int igusa_I6prime_fmpz(fmpz_t I6prime, fmpz* I)
{
  fmpz_t res, temp;
  int r;

  fmpz_init(res);
  fmpz_init(temp);
  
  /* Get I6prime from I6 */
  fmpz_mul_si(res, &I[2], -3);
  fmpz_mul(temp, &I[0], &I[1]);
  fmpz_add(res, res, temp);
  r = fmpz_divisible_si(res, 2);
  if (r)
    {
      fmpz_divexact_si(res, res, 2);
      fmpz_set(I6prime, res);
    }
  
  fmpz_clear(res);
  fmpz_clear(temp);
  return r;
}
