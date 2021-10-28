
#include "siegel.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("sp2gz_is_correct....");
  fflush(stdout);

  flint_randinit(state);

  for (iter = 0; iter < 100 * arb_test_multiplier(); iter++)
    {
      sp2gz_t m;
      slong g;

      g = 1 + n_randint(state, 10);
      sp2gz_init(m, g);

      sp2gz_randtest_triangular(m, state, n_randint(state, 50));
      if (!sp2gz_is_correct(m))
        {
	  flint_printf("FAIL\n");
	  flint_printf("g = %wd\n", g);
	  flint_printf("m = "); sp2gz_print(m); flint_printf("\n");
	  flint_abort();
	}
      sp2gz_randtest_diagonal(m, state, n_randint(state, 50));
      if (!sp2gz_is_correct(m))
        {
	  flint_printf("FAIL\n");
	  flint_printf("g = %wd\n", g);
	  flint_printf("m = "); sp2gz_print(m); flint_printf("\n");
	  flint_abort();
	}
      sp2gz_randtest(m, state, n_randint(state, 50));
      if (!sp2gz_is_correct(m))
        {
	  flint_printf("FAIL\n");
	  flint_printf("g = %wd\n", g);
	  flint_printf("m = "); sp2gz_print(m); flint_printf("\n");
	  flint_abort();
	}

      sp2gz_clear(m);
    }
  
  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

