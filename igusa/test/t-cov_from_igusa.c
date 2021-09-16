
#include "igusa.h"

int main()
{
  
  slong iter;
  flint_rand_t state;
  
  flint_printf("cov_from_igusa....");
  fflush(stdout);

  flint_randinit(state);

  for (iter = 0; iter < 500 * arb_test_multiplier(); iter++)
    {
      acb_ptr j_test;
      acb_ptr I;
      acb_ptr j;
      slong k;
      slong prec = 10 + n_randint(state, 500);
      slong mag_bits = 1 + n_randint(state, 10);
      int res = 1;

      j_test = _acb_vec_init(3);
      I = _acb_vec_init(4);
      j = _acb_vec_init(3);

      for (k = 0; k < 3; k++)
	{
	  acb_randtest_precise(&j[k], state, prec, mag_bits);
	}
      cov_from_igusa(I, j, prec);
      igusa_from_cov(j_test, I, prec);

      for (k = 0; k < 3; k++)
	{
	  if (!acb_overlaps(&j[k], &j_test[k])) res = 0;
	}
      if (!res)
	{
	  flint_printf("FAIL (j != j_test)\n");
	  for (k = 0; k < 3; k++)
	    {
	      acb_printd(&j_test[k], 30);
	      flint_printf("\n");
	      acb_printd(&j[k], 30);
	      flint_printf("\n");
	      fflush(stdout);
	      flint_abort();
	    }
	}

      _acb_vec_clear(j_test, 3);
      _acb_vec_clear(I, 4);
      _acb_vec_clear(j, 3);
    }
  
  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

