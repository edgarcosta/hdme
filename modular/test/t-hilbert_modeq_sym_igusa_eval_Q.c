
#include "modular.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("siegel_modeq_eval_Q....");
  fflush(stdout);

  flint_randinit(state);

  for (iter = 0; iter < 1 * arb_test_multiplier(); iter++)
    {
      fmpz_poly_t num1, num2, num3;
      fmpz_t den;
      slong delta;
      slong ell;
      fmpz_poly_t beta;
      fmpq* rs;
      slong rs_bits = 1 + n_randint(state, 10);
      slong k;
      int res;
      slong delta_max = 100;
      slong ell_max = 50;

      fmpz_poly_init(num1);
      fmpz_poly_init(num2);
      fmpz_poly_init(num3);
      fmpz_init(den);
      fmpz_poly_init(beta);
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
			  fmpq_randtest_not_zero(&rs[k], state, rs_bits);
			}
		      flint_printf("delta = %wd; ell = %wd; parameters are\n", delta, ell);
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
      _fmpq_vec_clear(rs, 2);
    }

  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
