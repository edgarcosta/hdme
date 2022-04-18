
#include "igusa.h"

void cov_rescale_fmpz(fmpz* I, fmpz* S, const fmpz_t scal)
{
  fmpz* res;
  fmpz_t f;

  res = _fmpz_vec_init(4);
  fmpz_init(f);
  
  fmpz_set(f, scal);
  fmpz_mul(&res[0], &S[0], f);
  fmpz_set(f, scal);
  fmpz_pow_ui(f, f, 2);
  fmpz_mul(&res[1], &S[1], f);
  fmpz_set(f, scal);
  fmpz_pow_ui(f, f, 3);
  fmpz_mul(&res[2], &S[2], f);
  fmpz_set(f, scal);
  fmpz_pow_ui(f, f, 5);
  fmpz_mul(&res[3], &S[3], f);

  _fmpz_vec_set(I, res, 4);
  
  _fmpz_vec_clear(res, 4);
  fmpz_clear(f);
}
