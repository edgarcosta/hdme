
#include "hilbert.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("gundlach_from_hilbert_param....");
  fflush(stdout);

  flint_randinit(state);

  for (iter = 0; iter < 100 * arb_test_multiplier(); iter++)
    {
      slong delta = 5;
      slong prec = 1000 + n_randint(state, 2000);
      fmpq* mn;
      fmpq* g;
      slong mn_bits = 1 + n_randint(state, 10);
      acb_ptr rs;
      acb_ptr I;
      acb_ptr g_acb;
      acb_ptr j;
      acb_ptr j_test;
      slong k;
      int res = 1;

      mn = _fmpq_vec_init(2);
      g = _fmpq_vec_init(2);
      rs = _acb_vec_init(2);
      I = _acb_vec_init(4);
      g_acb = _acb_vec_init(2);
      j = _acb_vec_init(3);
      j_test = _acb_vec_init(3);

      for (k = 0; k < 2; k++) fmpq_randbits(&mn[k], state, mn_bits);
      gundlach_from_hilbert_param(g, mn, delta);
      for (k = 0; k < 2; k++) acb_set_fmpq(&g_acb[k], &g[k], prec);
      igusa_from_gundlach(j, g_acb, delta, prec);
      
      acb_set_fmpq(&rs[0], &mn[0], prec);
      acb_set_fmpq(&rs[1], &mn[1], prec);
      hilbert_parametrize(I, rs, delta, prec);
      igusa_from_cov(j_test, I, prec);

      for (k = 0; k < 3; k++)
	{
	  if (!acb_overlaps(&j[k], &j_test[k])) res = 0;
	}
      if (!res)
	{
	  flint_printf("FAIL\n");
	  for (k = 0; k < 3; k++)
	    {
	      acb_printd(&j[k], 30); flint_printf("\n");
	      acb_printd(&j_test[k], 30); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}
      
      _fmpq_vec_clear(mn, 2);
      _fmpq_vec_clear(g, 2);
      _acb_vec_clear(rs, 2);
      _acb_vec_clear(I, 4);
      _acb_vec_clear(g_acb, 2);
      _acb_vec_clear(j, 3);
      _acb_vec_clear(j_test, 3);

    }
  

  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}



      
      
