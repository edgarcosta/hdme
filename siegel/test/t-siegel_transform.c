
#include "siegel.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("siegel_transform....");
  fflush(stdout);
  
  flint_randinit(state);
  
  /* Check associativity */
  for (iter = 0; iter < 100 * arb_test_multiplier(); iter++)
    {
      sp2gz_t m, n, mn;
      acb_mat_t z1, z2, z3;
      slong g, prec;
      int valid_transform;

      g = 1 + n_randint(state, 10);
      sp2gz_init(m, g);
      sp2gz_init(n, g);
      sp2gz_init(mn, g);
      acb_mat_init(z1, g, g);
      acb_mat_init(z2, g, g);
      acb_mat_init(z3, g, g);

      sp2gz_randtest(m, state, n_randint(state, 20));
      sp2gz_randtest(n, state, n_randint(state, 20));
      sp2gz_mul(mn, m, n);

      prec = 200 + n_randint(state, 200);

      siegel_halfspace_randtest(z1, state, prec);
      valid_transform = siegel_transform(z2, mn, z1, prec);

      if (!valid_transform)
	{
	  flint_printf("FAIL (insufficient precision to invert)\n");
	  flint_printf("g = %wd\n", g);
	  flint_printf("mn = "); sp2gz_print(mn); flint_printf("\n\n");
	  flint_printf("z1 = "); acb_mat_printd(z1, 30); flint_printf("\n\n");
	  flint_abort();
	}

      valid_transform = siegel_transform(z3, n, z1, prec);
      
      if (!valid_transform)
	{
	  flint_printf("FAIL (insufficient precision to invert)\n");
	  flint_printf("g = %wd\n", g);
	  flint_printf("n = "); sp2gz_print(n); flint_printf("\n\n");
	  flint_printf("z1 = "); acb_mat_printd(z1, 30); flint_printf("\n\n");
	  flint_abort();
	}

      valid_transform = siegel_transform(z3, m, z3, prec);

      if (!valid_transform)
	{
	  flint_printf("FAIL (insufficient precision to invert)\n");
	  flint_printf("g = %wd\n", g);
	  flint_printf("m = "); sp2gz_print(m); flint_printf("\n\n");
	  flint_printf("z3 = "); acb_mat_printd(z3, 30); flint_printf("\n\n");
	  flint_abort();
	}

      if (!acb_mat_overlaps(z2, z3))
	{
	  
	  flint_printf("FAIL\n");
	  flint_printf("g = %wd\n", g);
	  flint_printf("m = "); sp2gz_print(m); flint_printf("\n\n");
	  flint_printf("n = "); sp2gz_print(n); flint_printf("\n\n");
	  flint_printf("z1 = "); acb_mat_printd(z1, 30); flint_printf("\n\n");
	  flint_printf("z2 = "); acb_mat_printd(z2, 30); flint_printf("\n\n");
	  flint_printf("z3 = "); acb_mat_printd(z3, 30); flint_printf("\n\n");
	  flint_abort();
	}

      sp2gz_clear(m);
      sp2gz_clear(n);
      sp2gz_clear(mn);
      acb_mat_clear(z1);
      acb_mat_clear(z2);
      acb_mat_clear(z3);
    }
  
  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
