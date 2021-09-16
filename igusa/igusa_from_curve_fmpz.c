
#include "igusa.h"

void igusa_from_curve_fmpz(fmpq* j, const fmpz_poly_t crv)
{
  fmpz* I;
  
  I = _fmpz_vec_init(4);

  igusa_scalar_covariants_fmpz(I, crv);
  igusa_from_cov_fmpz(j, I);
  
  _fmpz_vec_clear(I, 4);
}
