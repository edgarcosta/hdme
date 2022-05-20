
#include "igusa.h"

void cov_normalize_fmpz(fmpz* I, fmpz* S, slong nb, slong* weights)
{
  slong p, v;
  fmpz_t g, f;
  slong i0 = -1;
  slong j;
  slong max_p = COV_MAXP;
  
  fmpz_init(g);
  fmpz_init(f);
  _fmpz_vec_set(I, S, nb);
  
  /* Compute i0 st index is nonzero */
  for (j = 0; j < nb; j++)
    {
      if (!fmpz_is_zero(&I[j]))
	{
	  i0 = j;
	  break;
	}
    }
  if (i0 == -1)
    {
      flint_printf("(cov_normalize_fmpz) Error: all entries are zero\n");
      fflush(stdout);
      flint_abort();
    }
  
  /* Rescale using valuations at small primes */
  p = 2;
  while (p < max_p)
    {
      fmpz_set_si(f, p);
      v = fmpz_remove(g, &I[i0], f) / (weights[i0]);
      for (j = 0; j < nb; j++)
	{
	  if (j != i0 && !fmpz_is_zero(&I[j]))
	    {
	      v = FLINT_MIN(v, fmpz_remove(g, &I[j], f) / (weights[j]));
	    }
	}      
      fmpz_pow_ui(g, f, v);
      cov_divexact_fmpz(I, I, g, nb, weights);
      p = n_nextprime(p, 1);
    }
  
  fmpz_gcd(g, &I[0], &I[1]);
  fmpz_gcd(g, g, &I[2]);
  fmpz_gcd(g, g, &I[3]);

  /* Rescale using gcd if possible, ignoring small primes */
  p = 2;
  while (p < max_p)
    {
      fmpz_set_si(f, p);
      fmpz_remove(g, g, f);
      p = n_nextprime(p, 1);
    }
  
  if (cov_divisible_fmpz(I, g, nb, weights))
    {
      cov_divexact_fmpz(I, I, g, nb, weights);
    }
  
  /* Update a sign in case of odd weight */
  for (j = 0; j < nb; j++)
    {
      if (!fmpz_is_zero(&I[j]) && (weights[j] %2) == 1)
	{
	  if (fmpz_cmp_si(&I[j], 0) < 0) cov_rescale_fmpz_si(I, I, -1, nb, weights); 
	  break;
	}
    }

  fmpz_clear(f);
  fmpz_clear(g);
}
