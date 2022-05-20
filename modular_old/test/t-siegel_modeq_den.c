
#include "modular.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("siegel_modeq_den....");
  fflush(stdout);

  flint_randinit(state);

  for (iter = 0; iter < 1 * arb_test_multiplier(); iter++)
    {
      slong prec = 5000 + n_randint(state, 5000);
      slong j_bits = 4 + n_randint(state, 10);
      acb_ptr j;
      acb_ptr I;
      fmpz_t temp;
      acb_mat_t tau;
      slong ell = 2; /* n_randprime(state, 3, 1); */ /* 2, 3, 5 or 7 */
      slong n = siegel_nb_cosets(ell);
      acb_ptr th2_vec, I_vec;
      acb_ptr stardets;
      acb_ptr th2;
      acb_t scal;
      acb_t den;
      slong k;
      int res;

      j = _acb_vec_init(3);
      I = _acb_vec_init(4);
      fmpz_init(temp);
      acb_mat_init(tau, 2, 2);
      th2_vec = _acb_vec_init(16 * n);
      I_vec = _acb_vec_init(4 * n);
      stardets = _acb_vec_init(n);
      th2 = _acb_vec_init(16);
      acb_init(scal);
      acb_init(den);

      for (k = 0; k < 3; k++)
	{
	  fmpz_randtest_not_zero(temp, state, j_bits);
	  acb_set_fmpz(&j[k], temp);
	}
      cov_from_igusa(I, j, prec);
      
      res = tau_theta2_from_igusa(tau, th2, I, prec);
      if (!res)
	{
	  flint_printf("FAIL (tau)\n");
	  for (k = 0; k < 3; k++)
	    {
	      acb_printd(&j[k], 30); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}

      res = siegel_modeq_theta2(th2_vec, stardets, tau, ell, prec);
      if (!res)
	{
	  flint_printf("FAIL (theta)\n");
	  for (k = 0; k < 3; k++)
	    {
	      acb_printd(&j[k], 30); flint_printf("\n");
	    }
	  acb_mat_printd(tau, 10);
	  fflush(stdout);
	  flint_abort();
	}

      modeq_cov(I_vec, th2_vec, n, prec);
      cov_from_theta2(I, th2, prec);
      siegel_modeq_scalar(scal, I, stardets, ell, prec);
      
      siegel_modeq_den(den, I_vec, scal, ell, prec);
      res = modeq_round_coeff(temp, den);
      if (!res)
	{
	  flint_printf("FAIL (integer)\n");
	  acb_printd(den, 500);
	  flint_printf("\n");
	  
	  flint_printf("Theta[0..15]:\n");
	  for (k = 0; k < 16; k++)
	    {
	      acb_printd(&th2_vec[k], 30); flint_printf("\n");
	    }
	  flint_printf("I[0..3]:\n");
	  for (k = 0; k < 4; k++)
	    {
	      acb_printd(&I_vec[k], 30); flint_printf("\n");
	    }
	  flint_printf("I_tau:\n");
	  for (k = 0; k < 4; k++)
	    {
	      acb_printd(&I[k], 30); flint_printf("\n");
	    }
	  igusa_from_cov(j, I, prec);
	  flint_printf("j:\n");
	  for (k = 0; k < 3; k++)
	    {
	      acb_printd(&j[k], 30); flint_printf("\n");
	    }
	}
      
      _acb_vec_clear(j, 3);
      _acb_vec_clear(I, 4);
      fmpz_clear(temp);
      acb_mat_clear(tau);
      _acb_vec_clear(th2_vec, 16 * n);
      _acb_vec_clear(I_vec, 4 * n);
      _acb_vec_clear(stardets, n);
      _acb_vec_clear(th2, 16);
      acb_clear(scal);
      acb_clear(den);
    }
  
  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
      


      
