
#include "modular.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("hilbert_modeq_igusa_eval_Q....");
  fflush(stdout);

  flint_randinit(state);

  for (iter = 0; iter < 1 * arb_test_multiplier(); iter++)
    {
      fmpz_poly_struct num_vec[3];
      fmpz_t den;
      slong delta;
      slong ell;
      fmpz_poly_t beta;
      fmpz_poly_t betabar;
      fmpq* rs;
      slong rs_bits = 2 + n_randint(state, 5);
      slong k;
      int res;
      slong delta_max = 15;
      slong ell_max = 12;

      for (k = 0; k < 3; k++) fmpz_poly_init(&num_vec[k]);
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
		      res = hilbert_modeq_igusa_eval_Q(num_vec, den, rs, ell, delta);
		      /* fmpz_print(den); flint_printf("\n");*/
		      /* fmpz_poly_print_pretty(num1, "x"); flint_printf("\n"); */
		      if (!res)
			{
			  flint_printf("FAIL\n");
			  fflush(stdout);
			  flint_abort();
			}
		      flint_printf("Denominator is a %wd-bit integer\n", 
				   fmpz_bits(den));
		    }
		}
	    }
	}
      
      for (k = 0; k < 3; k++) fmpz_poly_clear(&num_vec[k]);
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
