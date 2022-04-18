
#include "igusa.h"

void cov_normalize_fmpz(fmpz* I, fmpz* S)
{
  slong p, v;
  fmpz_t g, f;
  
  fmpz_init(g);
  fmpz_init(f);

  _fmpz_vec_set(I, S, 4);
  /* Rescale using valuations at small primes, then gcd */
  for (p = 2; p < 10000; p++)
    {
      if (n_is_prime(p))
	{
	  fmpz_set_si(f, p);
	  v = fmpz_remove(g, &I[1], f)/2;
	  if (!fmpz_is_zero(&I[0])) v = FLINT_MIN(v, fmpz_remove(g, &I[0], f));
	  if (!fmpz_is_zero(&I[2])) v = FLINT_MIN(v, fmpz_remove(g, &I[2], f)/3);
	  v = FLINT_MIN(v, fmpz_remove(g, &I[3], f)/5);
	  /* if (v != 0) flint_printf("Valuation of %d: %d\n", p, v); */

	  fmpz_pow_ui(g, f, v);
	  cov_divexact_fmpz(I, I, g);
	}
    }
  
  fmpz_gcd(g, &I[0], &I[1]);
  for (p = 2; p < 10000; p++)
    {
      fmpz_set_si(f, p);
      fmpz_remove(g, g, f);
    }
  if (cov_divisible_fmpz(I, g)) cov_divexact_fmpz(I, I, g);
  if (fmpz_cmp_si(&I[0], 0) < 0) cov_rescale_fmpz_si(I, I, -1);  

  fmpz_clear(f);
  fmpz_clear(g);
}
