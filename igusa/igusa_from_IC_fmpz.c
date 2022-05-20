
#include "igusa.h"

void igusa_from_IC_fmpz(fmpz* I, fmpz* IC)
{
  fmpz_t I6prime, I12;
  fmpz* resc;
  int r;

  fmpz_init(I6prime);
  fmpz_init(I12);
  resc = _fmpz_vec_init(4);
  _fmpz_vec_set(resc, IC, 4);
  
  r = igusa_I6prime_fmpz(I6prime, IC);
  if (!r)
    {
      fmpz_mul_si(&resc[0], &IC[0], n_pow(2, 1));
      fmpz_mul_si(&resc[1], &IC[1], n_pow(2, 2));
      fmpz_mul_si(&resc[2], &IC[2], n_pow(2, 3));
      fmpz_mul_si(&resc[3], &IC[3], n_pow(2, 5));      
    }
  igusa_I6prime_fmpz(I6prime, resc);
  fmpz_mul(I12, &resc[0], &resc[3]);

  fmpz_set(cov_I4(I), &resc[1]);
  fmpz_set(cov_I6prime(I), I6prime);
  fmpz_set(cov_I10(I), &resc[3]);
  fmpz_set(cov_I12(I), I12);

  cov_normalize_fmpz(I, I);

  fmpz_clear(I6prime);
  fmpz_clear(I12);
  _fmpz_vec_clear(resc, 4);  
}
