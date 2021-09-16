
#include "igusa.h"

void igusa_scalar_covariants_fmpz(fmpz* I, const fmpz_poly_t crv)
{
  /* Support weird aliasing */
  fmpz* res;
  fmpz_t a0, a1, a2, a3, a4, a5, a6;
  res = _fmpz_vec_init(4);

  fmpz_init(a0);
  fmpz_init(a1);
  fmpz_init(a2);
  fmpz_init(a3);
  fmpz_init(a4);
  fmpz_init(a5);
  fmpz_init(a6);

  fmpz_poly_get_coeff_fmpz(a0, crv, 0);
  fmpz_poly_get_coeff_fmpz(a1, crv, 1);
  fmpz_poly_get_coeff_fmpz(a2, crv, 2);
  fmpz_poly_get_coeff_fmpz(a3, crv, 3);
  fmpz_poly_get_coeff_fmpz(a4, crv, 4);
  fmpz_poly_get_coeff_fmpz(a5, crv, 5);
  fmpz_poly_get_coeff_fmpz(a6, crv, 6);

  igusa_I2_fmpz(&res[0], a0, a1, a2, a3, a4, a5, a6);
  igusa_I4_fmpz(&res[1], a0, a1, a2, a3, a4, a5, a6);
  igusa_I6prime_fmpz(&res[2], a0, a1, a2, a3, a4, a5, a6);
  igusa_I10_fmpz(&res[3], a0, a1, a2, a3, a4, a5, a6);
  _fmpz_vec_set(I, res, 4);
  
  _fmpz_vec_clear(res, 4);
  fmpz_clear(a0);
  fmpz_clear(a1);
  fmpz_clear(a2);
  fmpz_clear(a3);
  fmpz_clear(a4);
  fmpz_clear(a5);
  fmpz_clear(a6);
}
