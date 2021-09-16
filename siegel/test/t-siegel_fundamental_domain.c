
#include "siegel.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("siegel_fundamental_domain....");
  fflush(stdout);
  
  flint_randinit(state);

  for (iter = 0; iter < 100 * arb_test_multiplier(); iter++)
    {
      slong g = 2;
      slong bits = 10 + n_randint(state, 500);
      slong prec = 10 * g * bits;
      int res;

      arb_t tol;
      acb_mat_t z;
      acb_mat_t w;
      sp2gz_t m;

      arb_init(tol);
      acb_mat_init(z, g, g);
      acb_mat_init(w, g, g);
      sp2gz_init(m, g);
      
      arb_set_si(tol, 1);
      arb_mul_2exp_si(tol, tol, -bits); /* tol is larger than 2^(-prec) */

      siegel_halfspace_randtest(z, state, prec);
      /* z is not too far from the fundamental domain. */
      res = siegel_fundamental_domain(z, m, z, tol, prec);
      
      if (!res || !siegel_is_in_fundamental_domain(z, tol, prec)
	  || !sp2gz_is_correct(m))
	{ 
	  flint_printf("FAIL\n");
	  flint_printf("prec = %wd\n", prec);
	  flint_printf("z = "); acb_mat_printd(z, 30); flint_printf("\n\n");
	  flint_printf("m = "); sp2gz_print(m); flint_printf("\n\n");
	  flint_abort();
	}

      /* Random action of Sp_{2g}(ZZ) */
      sp2gz_randtest(m, state, bits);
      res = siegel_transform(w, m, z, prec);
      sp2gz_inv(m, m);
      res = res && siegel_transform(z, m, w, prec);
      
      if (!res)
	{
	  flint_printf("FAIL\n");
	  flint_printf("prec = %wd\n", prec);
	  flint_printf("z = "); acb_mat_printd(z, 30); flint_printf("\n\n");
	  flint_printf("m = "); sp2gz_print(m); flint_printf("\n\n");
	  flint_printf("w = "); acb_mat_printd(w, 30); flint_printf("\n\n");
	  flint_abort();
	}

      res = siegel_fundamental_domain(z, m, w, tol, prec);
      
      if (!res || !siegel_is_in_fundamental_domain(z, tol, prec)
	  || !sp2gz_is_correct(m))
	{ 
	  flint_printf("FAIL\n");
	  flint_printf("res = %wd\n", res);
	  flint_printf("prec = %wd\n", prec);
	  flint_printf("w = "); acb_mat_printd(w, 30); flint_printf("\n\n");
	  flint_printf("m = "); sp2gz_print(m); flint_printf("\n\n");
	  flint_printf("z = "); acb_mat_printd(z, 30); flint_printf("\n\n");
	  flint_abort();
	}
      
      arb_clear(tol);
      acb_mat_clear(z);
      acb_mat_clear(w);
      sp2gz_clear(m);
    }

  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
