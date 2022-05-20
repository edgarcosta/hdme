
#include "igusa.h"

void igusa_from_IC_fmpz(fmpz* I, fmpz* IC)
{
  fmpz_t I6prime, I12;
  slong weights[4] = IGUSA_HALFWEIGHTS;
  fmpz* resc;

  fmpz_init(I6prime);
  fmpz_init(I12);
  resc = _fmpz_vec_init(4);
  _fmpz_vec_set(resc, IC, 4);

  cov_rescale_fmpz_si(resc, IC, 2, 4, weights);
  igusa_I6prime_from_IC_fmpz(I6prime, resc);
  fmpz_mul(I12, &resc[0], &resc[3]);

  fmpz_set(igusa_I4(I), &resc[1]);
  fmpz_set(igusa_I6prime(I), I6prime);
  fmpz_set(igusa_I10(I), &resc[3]);
  fmpz_set(igusa_I12(I), I12);

  cov_normalize_fmpz(I, I, 4, weights);

  fmpz_clear(I6prime);
  fmpz_clear(I12);
  _fmpz_vec_clear(resc, 4);  
}
