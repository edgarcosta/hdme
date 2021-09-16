
#include "modular.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("siegel_modeq_eval_Fp....");
  fflush(stdout);

  flint_randinit(state);

  for (iter = 0; iter < 1 * arb_test_multiplier(); iter++)
    {
      fmpz_mod_poly_t pol1, pol2, pol3;
      fmpz_mod_ctx_t ctx;
      fmpz_t p;
      fmpz* j;
      slong ell =  2; /* n_randprime(state, 3, 1); */ /* 2, 3, 5 or 7 */
      slong p_bits = 1 + n_randint(state, 100);
      slong k;
      int res;

      fmpz_init(p);
      fmpz_randprime(p, state, p_bits, 1);
      fmpz_mod_ctx_init(ctx, p);
      
      fmpz_mod_poly_init(pol1, ctx);
      fmpz_mod_poly_init(pol2, ctx);
      fmpz_mod_poly_init(pol3, ctx);
      j = _fmpz_vec_init(3);

      for (k = 0; k < 3; k++)
	{
	  do
	    {
	      fmpz_randtest_mod(&j[k], state, p);
	    }
	  while (fmpz_divisible(&j[k], p)); /* Nonzero mod p */
	}
      
      res = siegel_modeq_eval_Fp(pol1, pol2, pol3, j, ell, ctx);
      if (!res)
	{
	  flint_printf("FAIL\n");
	  for (k = 0; k < 3; k++)
	    {
	      flint_printf("p = ");
	      fmpz_print(p); flint_printf("\n");
	      fmpz_print(&j[k]); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}
      
      fmpz_mod_poly_clear(pol1, ctx);
      fmpz_mod_poly_clear(pol2, ctx);
      fmpz_mod_poly_clear(pol3, ctx);
      fmpz_clear(p);
      fmpz_mod_ctx_clear(ctx);
      _fmpz_vec_clear(j, 3);
    }
  
  
  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
