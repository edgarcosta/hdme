
#include "modular.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("siegel_modeq_num....");
  fflush(stdout);

  flint_randinit(state);

  for (iter = 0; iter < 1 * arb_test_multiplier(); iter++)
    {
      slong prec = 10000 + n_randint(state, 5000);
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
      acb_poly_t num1_acb, num2_acb, num3_acb;
      acb_t den_acb;
      fmpz_poly_t num1, num2, num3;
      fmpz_t den;
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
      acb_poly_init(num1_acb);
      acb_poly_init(num2_acb);
      acb_poly_init(num3_acb);
      acb_init(den_acb);
      fmpz_poly_init(num1);
      fmpz_poly_init(num2);
      fmpz_poly_init(num3);
      fmpz_init(den);

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

      siegel_modeq_cov(I_vec, th2_vec, ell, prec);
      cov_from_theta2(I, th2, prec);
      siegel_modeq_scalar(scal, I, stardets, ell, prec);
      
      siegel_modeq_den(den_acb, I_vec, scal, ell, prec);
      siegel_modeq_num(num1_acb, num2_acb, num3_acb, I_vec, scal, ell, prec);
      
      res = siegel_modeq_round(num1, num2, num3, den,
			       num1_acb, num2_acb, num3_acb, den_acb, ell);
      if (!res)
	{
	  flint_printf("FAIL (integers)\n");
	  acb_poly_printd(num1_acb, 500);
	}
      res = ((fmpz_poly_degree(num1) == n)
	     && (fmpz_poly_degree(num2) == n-1)
	     && (fmpz_poly_degree(num3) == n-1));
      if (!res)
	{
	  flint_printf("FAIL (degrees)\n");
	  fmpz_poly_print_pretty(num1, "X");
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
      acb_poly_clear(num1_acb);
      acb_poly_clear(num2_acb);
      acb_poly_clear(num3_acb);
      acb_clear(den_acb);
      fmpz_poly_clear(num1);
      fmpz_poly_clear(num2);
      fmpz_poly_clear(num3);
      fmpz_clear(den);
    }
  
  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
      


      
