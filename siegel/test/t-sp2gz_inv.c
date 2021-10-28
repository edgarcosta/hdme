
#include "siegel.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("sp2gz_inv....");
  fflush(stdout);
  
  flint_randinit(state);

  for (iter = 0; iter < 100 * arb_test_multiplier(); iter++)
    {
      sp2gz_t m, n, p, i;
      slong g = 1 + n_randint(state, 10);

      sp2gz_init(m, g);
      sp2gz_init(n, g);
      sp2gz_init(p, g);
      sp2gz_init(i, g);

      sp2gz_randtest(m, state, n_randint(state, 50));
      sp2gz_randtest(n, state, n_randint(state, 50));
      sp2gz_randtest(p, state, n_randint(state, 50));

      sp2gz_inv(n, m);
      sp2gz_mul(p, m, n);
      sp2gz_one(i);

      if (!sp2gz_equal(p, i) || !sp2gz_is_correct(n))
        {
	  flint_printf("FAIL\n");
	  flint_printf("m = "); sp2gz_print(m); flint_printf("\n");
	  flint_printf("n = "); sp2gz_print(n); flint_printf("\n");
	  flint_printf("p = "); sp2gz_print(p); flint_printf("\n");
	  flint_printf("i = "); sp2gz_print(i); flint_printf("\n");
	  flint_abort();
        }

      sp2gz_inv(m, m);

      if (!sp2gz_equal(m, n) || !sp2gz_is_correct(m))
        {
	  flint_printf("FAIL (aliasing)\n");
	  flint_printf("m = "); sp2gz_print(m); flint_printf("\n");
	  flint_printf("n = "); sp2gz_print(n); flint_printf("\n");
	  flint_abort();
        }

      sp2gz_clear(m);
      sp2gz_clear(n);
      sp2gz_clear(p);
      sp2gz_clear(i);
    }

  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
