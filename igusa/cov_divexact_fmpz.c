
#include "igusa.h"

void cov_divexact_fmpz(fmpz* I, fmpz* S, const fmpz_t scal)
{
  fmpz* res;
  fmpz_t f;

  res = _fmpz_vec_init(4);
  fmpz_init(f);

  fmpz_divexact(&res[0], &S[0], scal);
  fmpz_pow_ui(f, scal, 2);
  fmpz_divexact(&res[1], &S[1], f);
  fmpz_pow_ui(f, scal, 3);
  fmpz_divexact(&res[2], &S[2], f);
  fmpz_pow_ui(f, scal, 5);
  fmpz_divexact(&res[3], &S[3], f);

  _fmpz_vec_set(I, res, 4);
  
  _fmpz_vec_clear(res, 4);
  fmpz_clear(f);
}
