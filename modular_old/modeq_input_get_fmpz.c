
#include "modular.h"

void modeq_input_get_fmpz(fmpz_t den, fmpz* num, fmpq* j, slong len)
{
  mpz_t n, d;
  fmpz_t aux, res;
  slong k;

  mpz_init(n);
  mpz_init(d);
  fmpz_init(aux);
  fmpz_init(res);

  fmpz_one(res);
  for (k = 0; k < len; k++)
    {
      fmpq_get_mpz_frac(n, d, &j[k]);
      fmpz_set_mpz(aux, d);
      fmpz_lcm(res, res, aux);
    }
  fmpz_set(den, res);
  
  for (k = 0; k < len; k++)
    {
      fmpq_get_mpz_frac(n, d, &j[k]);
      fmpz_set_mpz(aux, d);
      fmpz_divexact(res, den, aux);
      fmpz_set_mpz(aux, n);
      fmpz_mul(&num[k], aux, res);
    }
  
  fmpz_clear(aux);
  fmpz_clear(res);
  mpz_clear(n);
  mpz_clear(d);
}
