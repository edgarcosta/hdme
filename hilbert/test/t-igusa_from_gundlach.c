
#include "hilbert.h"


int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("igusa_from_gundlach....");
  fflush(stdout);

  flint_randinit(state);

  for (iter = 0; iter < 500 * arb_test_multiplier(); iter++)
    {
      slong delta = 5;
      acb_ptr j, g, I, test;
      acb_t r, s;
      slong prec = 100 + n_randint(state, 500);
      slong mag_bits = 1 + n_randint(state, 5);
      slong k;
      int res = 1;

      j = _acb_vec_init(3);
      g = _acb_vec_init(2);
      I = _acb_vec_init(4);
      test = _acb_vec_init(4);
      acb_init(r);
      acb_init(s);

      /* I (from humbert_parametrize) -> j -> g -> j */
      acb_randtest_precise(r, state, prec, mag_bits);
      acb_randtest_precise(s, state, prec, mag_bits);
      humbert_parametrize(I, r, s, delta, prec);
      igusa_from_cov(j, I, prec);
      gundlach_from_igusa(g, I, delta, prec);
      igusa_from_gundlach(test, g, delta, prec);
      
      for (k = 0; k < 3; k++)
	{
	  if (!acb_overlaps(&j[k], &test[k])) res = 0;
	}
      if (!res)
	{
	  flint_printf("FAIL (Igusa)\n");
	  fflush(stdout);
	  flint_abort();
	}

      /* g -> j -> g */
      acb_randtest_precise(&g[0], state, prec, mag_bits);
      acb_randtest_precise(&g[1], state, prec, mag_bits);
      igusa_from_gundlach(j, g, delta, prec);
      gundlach_from_igusa(test, j, delta, prec);
      
      for (k = 0; k < 2; k++)
	{
	  if (!acb_overlaps(&g[k], &test[k])) res = 0;
	}
      if (!res)
	{
	  flint_printf("FAIL (Gundlach)\n");
	  fflush(stdout);
	  flint_abort();
	}      

      _acb_vec_clear(j, 3);
      _acb_vec_clear(g, 2);
      _acb_vec_clear(I, 4);
      _acb_vec_clear(test, 4);
      acb_clear(r);
      acb_clear(s);      
    }
  
  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}


