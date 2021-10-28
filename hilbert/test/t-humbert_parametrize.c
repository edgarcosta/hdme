
#include "hilbert.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("humbert_parametrize....");
  fflush(stdout);

  flint_randinit(state);

  for (iter = 0; iter < 1 * arb_test_multiplier(); iter++)
    {
      slong delta;
      fmpq* rs;
      acb_t r_acb, s_acb;
      acb_ptr I, j, j_test;
      acb_t t1, t2;
      acb_mat_t tau;
      fmpz_mat_t eta;
      slong rs_bits = 5 + n_randint(state, 10);
      slong prec = 1000;
      slong k;
      int res;
      slong delta_max = 100;
      int v = 0;
      
      rs = _fmpq_vec_init(2);
      acb_init(r_acb);
      acb_init(s_acb);
      I = _acb_vec_init(4);
      j = _acb_vec_init(3);
      j_test = _acb_vec_init(3);
      acb_init(t1);
      acb_init(t2);
      acb_mat_init(tau, 2, 2);
      fmpz_mat_init(eta, 4, 4);
      
      for (delta = 5; delta < delta_max; delta++)
	{
	  if (hilbert_is_fundamental(delta))
	    {
	      for (k = 0; k < 2; k++)
		{
		  fmpq_randbits(&rs[k], state, rs_bits);
		}
	      if (v)
		{
		  flint_printf("delta = %wd; parameters are\n", delta);
		  fmpq_print(&rs[0]); flint_printf("\n");
		  fmpq_print(&rs[1]); flint_printf("\n");
		}
	      acb_set_fmpq(r_acb, &rs[0], prec);
	      acb_set_fmpq(s_acb, &rs[1], prec);
	      humbert_parametrize(I, r_acb, s_acb, delta, prec);
	      res = tau_from_igusa(tau, I, prec);
	      igusa_from_cov(j, I, prec);

	      if (!res)
		{
		  flint_printf("FAIL (period matrix)\n");
		  for (k = 0; k < 2; k++)
		    {
		      fmpq_print(&rs[k]); flint_printf("\n");
		    }
		  for (k = 0; k < 4; k++)
		    {
		      acb_printd(&I[k], 30); flint_printf("\n");
		    }
		  fflush(stdout);
		  flint_abort();
		}

	      res = hilbert_inverse(t1, t2, eta, tau, delta, prec);
	      
	      if (!res)
		{
		  flint_printf("FAIL (Hilbert inversion)\n");
		  for (k = 0; k < 2; k++)
		    {
		      fmpq_print(&rs[k]); flint_printf("\n");
		    }
		  fflush(stdout);
		  flint_abort();
		}

	      hilbert_map(tau, t1, t2, delta, prec);
	      igusa_from_tau(j_test, tau, prec);

	      for (k = 0; k < 2; k++)
		{
		  if (!acb_overlaps(&j_test[k], &j[k])) res = 0;
		}
	      
	      if (!res)
		{
		  flint_printf("FAIL (overlap)\n");
		  for (k = 0; k < 2; k++)
		    {
		      fmpq_print(&rs[k]); flint_printf("\n");
		    }
		  fflush(stdout);
		  flint_abort();
		}
	    }
	}

      _fmpq_vec_clear(rs, 2);
      acb_clear(r_acb);
      acb_clear(s_acb);
      _acb_vec_clear(I, 4);
      _acb_vec_clear(j, 3);
      _acb_vec_clear(j_test, 3);
      fmpz_mat_clear(eta);
    }

  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
