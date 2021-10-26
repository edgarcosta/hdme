
#include "modular.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("hilbert_modeq_gundlach_eval_Q....");
  fflush(stdout);

  flint_randinit(state);

  for (iter = 0; iter < 1 * arb_test_multiplier(); iter++)
    {
      fmpz_poly_t num1, num2;
      fmpz_t den;
      slong ell;
      fmpz_poly_t beta;
      fmpz_poly_t betabar;
      fmpq* g;
      slong g_bits;
      fmpz_t gnum, gden;
      slong k;
      int res;
      slong delta = 5;
      slong ell_max = 500;

      fmpz_poly_init(num1);
      fmpz_poly_init(num2);
      fmpz_init(den);
      fmpz_poly_init(beta);
      fmpz_poly_init(betabar);
      g = _fmpq_vec_init(2);
      fmpz_init(gnum);
      fmpz_init(gden);

      for (ell = 2; ell < ell_max; ell++)
	{
	  if (n_is_prime(ell) && hilbert_splits(beta, ell, delta))
	    {
	      g_bits = 5 + n_randint(state, 5);
	      for (k = 0; k < 2; k++)
		{
		  fmpz_randbits(gnum, state, g_bits);
		  fmpz_one(gden);
		  fmpq_set_fmpz_frac(&g[k], gnum, gden);
		}
	      hilbert_conjugate(betabar, beta, delta);
	      flint_printf("delta = %wd; ell = %wd; beta = ", delta, ell);
	      fmpz_poly_print_pretty(beta, "x");
	      flint_printf(", betabar = ");
	      fmpz_poly_print_pretty(betabar, "x");
	      flint_printf(", parameters are\n");		      
	      fmpq_print(&g[0]); flint_printf("\n");
	      fmpq_print(&g[1]); flint_printf("\n");		      
	      res = hilbert_modeq_gundlach_eval_Q(num1, num2, den, g, ell, delta);
	      /* fmpz_print(den); flint_printf("\n");*/
	      /* fmpz_poly_print_pretty(num1, "x"); flint_printf("\n"); */
	      if (!res)
		{
		  flint_printf("FAIL\n");
		  fflush(stdout);
		  flint_abort();
		}
	      flint_printf("Denominator is a %wd-bit integer\n", 
			   fmpz_bits(den));
	    }
	}

      fmpz_poly_clear(num1);
      fmpz_poly_clear(num2);
      fmpz_clear(den);
      fmpz_poly_clear(beta);
      fmpz_poly_clear(betabar);
      _fmpq_vec_clear(g, 2);
      fmpz_clear(gnum);
      fmpz_clear(gden);
    }

  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}


