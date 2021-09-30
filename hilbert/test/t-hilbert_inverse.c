
#include "hilbert.h"

int main()
{
  slong iter;
  flint_rand_t state;
  slong delta;
  
  flint_printf("hilbert_inverse....");
  fflush(stdout);
  
  flint_randinit(state);
  
  for (delta = 5; delta < 100; delta++)
    {
      if (hilbert_is_fundamental(delta))
	{
	  for (iter = 0; iter < 10 * arb_test_multiplier(); iter++)
	    {
	      acb_t t1, t2, t1_test, t2_test;
	      acb_mat_t tau, tau_test;
	      slong prec = 1000 + n_randint(state, 1000);
	      slong m_bits = 4;
	      sp2gz_t m;
	      sp2gz_t minv;
	      int res;
	      
	      acb_init(t1);
	      acb_init(t2);
	      acb_init(t1_test);
	      acb_init(t2_test);
	      acb_mat_init(tau, 2, 2);
	      acb_mat_init(tau_test, 2, 2);
	      sp2gz_init(m, 2);
	      sp2gz_init(minv, 2);
	      
	      hilbert_halfspace_randtest(t1, t2, state, prec);
	      hilbert_map(tau, t1, t2, delta, prec);
	      res = hilbert_inverse(t1_test, t2_test, minv, tau, delta, prec);

	      sp2gz_inv(m, minv);
	      hilbert_map(tau_test, t1_test, t2_test, delta, prec);
	      siegel_transform(tau_test, m, tau_test, prec);

	      if (!res || !acb_mat_overlaps(tau, tau_test))
		{
		  flint_printf("FAIL (wrong preimage)\n");
		  flint_printf("delta = %wd\n", delta);
		  flint_printf("(t1, t2):\n");
		  acb_printd(t1, 30); flint_printf("\n");
		  acb_printd(t2, 30); flint_printf("\n");
		  acb_mat_printd(tau, 10); flint_printf("\n");
		  flint_printf("Preimage:\n");
		  acb_printd(t1_test, 30); flint_printf("\n");
		  acb_printd(t2_test, 30); flint_printf("\n");
		  sp2gz_print(m);
		  acb_mat_printd(tau_test, 10); flint_printf("\n");
		  fflush(stdout);
		  flint_abort();
		}

	      sp2gz_randtest(m, state, m_bits);
	      siegel_transform(tau, m, tau, prec);
	      res = hilbert_inverse(t1_test, t2_test, minv, tau, delta, prec);
	      
	      sp2gz_inv(minv, minv);
	      hilbert_map(tau_test, t1_test, t2_test, delta, prec);
	      siegel_transform(tau_test, minv, tau_test, prec);
	      
	      if (!res || !acb_mat_overlaps(tau, tau_test))
		{
		  flint_printf("FAIL (wrong preimage)\n");
		  flint_printf("delta = %wd\n", delta);
		  flint_printf("(t1, t2):\n");
		  acb_printd(t1, 30); flint_printf("\n");
		  acb_printd(t2, 30); flint_printf("\n");
		  sp2gz_print(m);
		  acb_mat_printd(tau, 10); flint_printf("\n");
		  flint_printf("Preimage:\n");
		  acb_printd(t1_test, 30); flint_printf("\n");
		  acb_printd(t2_test, 30); flint_printf("\n");
		  sp2gz_print(minv);
		  acb_mat_printd(tau_test, 10); flint_printf("\n");
		  fflush(stdout);
		  flint_abort();
		}
	      
	      acb_clear(t1);
	      acb_clear(t2);
	      acb_clear(t1_test);
	      acb_clear(t2_test);
	      acb_mat_clear(tau);
	      acb_mat_clear(tau_test);
	      sp2gz_clear(m);
	      sp2gz_clear(minv);
	    }
	}
    }

  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

