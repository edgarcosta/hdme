
#include "hilbert.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("gundlach_cov_from_igusa....");
  fflush(stdout);

  flint_randinit(state);

  for (iter = 0; iter < 500 * arb_test_multiplier(); iter++)
    {
      slong delta = 5;
      acb_ptr j, g, G, I, test;
      acb_t r, s;
      slong prec = 100 + n_randint(state, 500);
      slong mag_bits = 1 + n_randint(state, 5);      
      slong k;
      int res = 1;

      j = _acb_vec_init(3);
      g = _acb_vec_init(2);
      G = _acb_vec_init(3);
      I = _acb_vec_init(4);
      test = _acb_vec_init(4);
      acb_init(r);
      acb_init(s);

      /* I (from humbert_parametrize) -> j -> g -> j */
      acb_randtest_precise(r, state, prec, mag_bits);
      acb_randtest_precise(s, state, prec, mag_bits);
      humbert_parametrize(I, r, s, delta, prec);
      
      gundlach_cov_from_igusa(G, I, delta, prec);
      gundlach_from_cov(g, G, delta, prec);
      gundlach_from_igusa(test, I, delta, prec);
      
      for (k = 0; k < 2; k++)
	{
	  if (!acb_overlaps(&g[k], &test[k])) res = 0;
	}
      if (!res)
	{
	  flint_printf("FAIL (Gundlach)\n");
	  flint_printf("Igusa covariants:\n");
	  for (k = 0; k < 4; k++)
	    {
	      acb_printd(&I[k], 30); flint_printf("\n");
	    }
	  flint_printf("Gundlach covariants:\n");
	  for (k = 0; k < 3; k++)
	    {
	      acb_printd(&G[k], 30); flint_printf("\n");
	    }	  
	  flint_printf("Gundlach invariants:\n");
	  for (k = 0; k < 2; k++)
	    {
	      acb_printd(&g[k], 30); flint_printf("\n");
	    }	  
	  flint_printf("Igusa invariants:\n");
	  for (k = 0; k < 3; k++)
	    {
	      acb_printd(&j[k], 30); flint_printf("\n");
	    }
	  flint_printf("Test (g-invariants):\n");
	  for (k = 0; k < 2; k++)
	    {
	      acb_printd(&test[k], 30); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}      

      _acb_vec_clear(j, 3);
      _acb_vec_clear(g, 2);
      _acb_vec_clear(I, 4);
      _acb_vec_clear(G, 3);
      _acb_vec_clear(test, 4);
      acb_clear(r);
      acb_clear(s);      
    }
  
  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
