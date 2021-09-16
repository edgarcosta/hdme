
#include "siegel.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("sp2gz_mul....");
  fflush(stdout);

  flint_randinit(state);

  for (iter = 0; iter < 500 * arb_test_multiplier(); iter++)
    {
      sp2gz_t m, n, p, u, v;
      slong g;

      g = 1 + n_randint(state, 10);
      sp2gz_init(m, g);
      sp2gz_init(n, g);
      sp2gz_init(p, g);
      sp2gz_init(u, g);
      sp2gz_init(v, g);
      
      sp2gz_randtest(m, state, n_randint(state, 50));
      sp2gz_randtest(n, state, n_randint(state, 50));
      sp2gz_randtest(p, state, n_randint(state, 50));
      sp2gz_randtest(u, state, n_randint(state, 50));
      sp2gz_randtest(v, state, n_randint(state, 50));
      
      /* test (m*n)*p = m*(n*p) */
      
      sp2gz_mul(u, m, n);
      sp2gz_mul(u, u, p);
      
      sp2gz_mul(v, n, p);
      sp2gz_mul(v, m, v);
      
      if (!sp2gz_equal(u, v) || !sp2gz_is_correct(u) || !sp2gz_is_correct(v))
        {
	  flint_printf("FAIL\n");
	  flint_printf("g = %wd\n", g);
	  flint_printf("m = "); sp2gz_print(m); flint_printf("\n");
	  flint_printf("n = "); sp2gz_print(n); flint_printf("\n");
	  flint_printf("p = "); sp2gz_print(p); flint_printf("\n");
	  flint_printf("u = "); sp2gz_print(u); flint_printf("\n");
	  flint_printf("v = "); sp2gz_print(v); flint_printf("\n");
	  flint_abort();
        }
      
      sp2gz_clear(m);
      sp2gz_clear(n);
      sp2gz_clear(p);
      sp2gz_clear(u);
      sp2gz_clear(v);
    }
  
  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
