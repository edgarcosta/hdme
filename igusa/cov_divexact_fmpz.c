
#include "igusa.h"

void cov_divexact_fmpz(fmpz* I, fmpz* S, const fmpz_t scal)
{
  fmpz* res;
  fmpz_t f;
  slong wt[4] = COV_WEIGHTS;
  slong j;

  res = _fmpz_vec_init(4);
  fmpz_init(f);

  for (j = 0; j < 4; j++)
    {
      fmpz_pow_ui(f, scal, wt[j]/2);
      fmpz_divexact(&res[j], &S[j], f);
    }

  _fmpz_vec_set(I, res, 4);
  
  _fmpz_vec_clear(res, 4);
  fmpz_clear(f);
}
