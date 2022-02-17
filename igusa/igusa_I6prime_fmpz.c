
#include "igusa.h"

void igusa_I6prime_fmpz(fmpz_t I6prime, fmpz* I)
{
  fmpz_t res, temp;

  fmpz_init(res);
  fmpz_init(temp);
  
  /* Get I6prime from I6 */
  fmpz_mul_si(res, &I[2], -3);
  fmpz_mul(temp, &I[0], &I[1]);
  fmpz_add(res, res, temp);
  fmpz_divexact_si(res, res, 2);

  fmpz_set(I6prime, res);
  fmpz_clear(res);
  fmpz_clear(temp);
}
