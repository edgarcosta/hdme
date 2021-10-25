
#include "modular.h"


int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("hilbert_modeq_gundlach_scalar....");
  fflush(stdout);

  flint_randinit(state);

  for (iter = 0; iter < 10 * arb_test_multiplier(); iter++)
    {
      acb_ptr g;
      fmpz_t g_fmpz;
      acb_ptr j;
      acb_ptr I;
      acb_ptr stardets;
      acb_t scal;
      acb_t temp;
      fmpz_t test;
      slong g_bits = 5 + n_randint(state, 20);
      slong delta = 5;
      slong ell = 2; /* Total weight is 2*10*3 = 60 */
      slong n = 2 * hilbert_nb_cosets(ell, delta);
      slong k;
      slong prec = 500 + n_randint(state, 1000);
      int res;

      g = _acb_vec_init(2);
      fmpz_init(g_fmpz);
      j = _acb_vec_init(3);
      I = _acb_vec_init(4);
      stardets = _acb_vec_init(n);
      acb_init(scal);
      acb_init(temp);
      fmpz_init(test);

      for (k = 0; k < 2; k++)
	{
	  fmpz_randbits(g_fmpz, state, g_bits);
	  acb_set_fmpz(&g[k], g_fmpz);
	}
      igusa_from_gundlach(j, g, delta, prec);
      cov_from_igusa(I, j, prec);

      for (k = 0; k < n; k++) acb_one(&stardets[k]);
      hilbert_modeq_gundlach_scalar(scal, I, stardets, ell, delta, prec);
      acb_pow_ui(temp, &I[2], 10, prec);
      acb_mul(scal, scal, temp, prec);
      res = siegel_modeq_round_coeff(test, scal);

      if (!res)
	{
	  flint_printf("FAIL\n");
	  acb_printd(scal, 30); flint_printf("\n");
	  flint_printf("Gundlach invariants:\n");
	  for (k = 0; k < 2; k++)
	    {
	      acb_printd(&g[k], 30);
	      flint_printf("\n");
	    }
	  flint_printf("Igusa covariants:\n");
	  for (k = 0; k < 4; k++)
	    {
	      acb_printd(&I[k], 30);
	      flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}
      
      _acb_vec_clear(g, 2);
      fmpz_clear(g_fmpz);
      _acb_vec_clear(j, 3);
      _acb_vec_clear(I, 4);
      _acb_vec_clear(stardets, n);
      acb_clear(scal);
      acb_clear(temp);
      fmpz_clear(test);
    }

  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
