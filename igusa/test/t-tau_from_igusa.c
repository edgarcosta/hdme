
#include "igusa.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("tau_from_igusa....");
  fflush(stdout);

  flint_randinit(state);

  for (iter = 0; iter < 10 * arb_test_multiplier(); iter++)
    {
      slong prec = 2000 + n_randint(state, 20000);
      acb_mat_t tau;
      acb_ptr j_test;
      acb_ptr j;
      acb_ptr I;
      slong mag_bits = 1 + n_randint(state, 10);
      slong k;
      int res;
      
      acb_mat_init(tau, 2, 2);
      j_test = _acb_vec_init(3);
      j = _acb_vec_init(3);
      I = _acb_vec_init(4);

      /* Generate Igusa invariants */
      for (k = 0; k < 3; k++) acb_randtest_precise(&j_test[k], state, prec, mag_bits);
      cov_from_igusa(I, j_test, prec);
      res = tau_from_igusa(tau, I, prec);

      if (!res)
	{ 
	  flint_printf("FAIL (tau)\n");
	  for (k = 0; k < 3; k++)
	    {
	      acb_printd(&j_test[k], 30); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}

      /* Compute Igusa invariants and check */
      res = igusa_from_tau(j, tau, prec);
      if (!res)
	{ 
	  flint_printf("FAIL (j)\n");
	  for (k = 0; k < 3; k++)
	    {
	      acb_printd(&j_test[k], 30); flint_printf("\n");
	    }
	  acb_mat_printd(tau, 10);
	  fflush(stdout);
	  flint_abort();
	}

      for (k = 0; k < 3; k++)
	{
	  if (!acb_overlaps(&j[k], &j_test[k])) res = 0;
	}
      if (!res)
	{
	  flint_printf("FAIL (overlap)\n");
	  for (k = 0; k < 3; k++)
	    {
	      acb_printd(&j_test[k], 30); flint_printf("\n");
	      acb_printd(&j[k], 30); flint_printf("\n");
	    }
	  acb_mat_printd(tau, 10);
	  fflush(stdout);
	  flint_abort();
	}

      acb_mat_clear(tau);
      _acb_vec_clear(j_test, 3);
      _acb_vec_clear(j, 3);
      _acb_vec_clear(I, 4);
    }
  
  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
      
