
#include "modular.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("hilbert_modeq_sym_igusa_eval_Q....");
  fflush(stdout);

  flint_randinit(state);

  for (iter = 0; iter < 1 * arb_test_multiplier(); iter++)
    {
      fmpz_poly_t num1, num2, num3;
      fmpz_t den;
      slong delta;
      slong ell;
      fmpz_poly_t beta;
      fmpz_poly_t betabar;
      fmpq* rs;
      slong rs_bits = 5 + n_randint(state, 10);
      slong k;
      int res;
      slong delta_max = 10;
      slong ell_max = 15;

      fmpz_poly_init(num1);
      fmpz_poly_init(num2);
      fmpz_poly_init(num3);
      fmpz_init(den);
      fmpz_poly_init(beta);
      fmpz_poly_init(betabar);
      rs = _fmpq_vec_init(2);
      
      for (delta = 5; delta < delta_max; delta++)
	{
	  if (hilbert_is_fundamental(delta))
	    {
	      for (ell = 2; ell < ell_max; ell++)
		{
		  if (n_is_prime(ell) && hilbert_splits(beta, ell, delta))
		    {
		      for (k = 0; k < 2; k++)
			{
			  fmpq_randbits(&rs[k], state, rs_bits);
			}
		      hilbert_conjugate(betabar, beta, delta);
		      flint_printf("delta = %wd; ell = %wd; beta = ", delta, ell);
		      fmpz_poly_print_pretty(beta, "x");
		      flint_printf(", betabar = ");
		      fmpz_poly_print_pretty(betabar, "x");
		      flint_printf(", parameters are\n");		      
		      fmpq_print(&rs[0]); flint_printf("\n");
		      fmpq_print(&rs[1]); flint_printf("\n");		      
		      res = hilbert_modeq_sym_igusa_eval_Q(num1, num2, num3, den, rs, ell, delta);
		      if (!res)
			{
			  flint_printf("FAIL\n");
			  for (k = 0; k < 2; k++)
			    {
			      fmpq_print(&rs[k]); flint_printf("\n");
			    }
			  fflush(stdout);
			  flint_abort();
			}
		    }
		}
	    }
	}
      
      fmpz_poly_clear(num1);
      fmpz_poly_clear(num2);
      fmpz_poly_clear(num3);
      fmpz_clear(den);
      fmpz_poly_clear(beta);
      fmpz_poly_clear(betabar);
      _fmpq_vec_clear(rs, 2);
    }

  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
