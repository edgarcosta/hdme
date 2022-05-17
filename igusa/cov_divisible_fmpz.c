
#include "igusa.h"

int cov_divisible_fmpz(fmpz* I, const fmpz_t scal)
{
  fmpz_t f;
  int r = 1;
  slong wt[4] = COV_WEIGHTS;
  slong j;

  fmpz_init(f);

  for (j = 0; j < 4; j++)
    {
      fmpz_pow_ui(f, scal, wt[j]/2);
      r = r && fmpz_divisible(&I[j], f);
    }

  fmpz_clear(f);
  return r;
}
