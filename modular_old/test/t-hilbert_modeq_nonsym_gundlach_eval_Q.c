
#include "modular.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("hilbert_modeq_nonsym_gundlach_eval_Q....");
  fflush(stdout);

  flint_randinit(state);

  for (iter = 0; iter < 1 * arb_test_multiplier(); iter++)
    {
      fmpz_poly_struct num_vec[2];
      fmpz_t den;
      slong ell;
      fmpz_poly_t beta;
      fmpq* mn;
      slong mn_bits;
      fmpz_t rnum, rden;
      slong k;
      int res;
      slong delta = 5;
      slong ell_max = 20;

      for (k = 0; k < 2; k++) fmpz_poly_init(&num_vec[k]);
      fmpz_init(den);
      fmpz_poly_init(beta);
      mn = _fmpq_vec_init(2);
      fmpz_init(rnum);
      fmpz_init(rden);
      
      for (ell = 2; ell < ell_max; ell++)
	{
	  if (n_is_prime(ell) && hilbert_splits(beta, ell, delta))
	    {
	      mn_bits = 2 + n_randint(state, 3);
	      for (k = 0; k < 2; k++)
		{
		  fmpz_randbits(rnum, state, mn_bits);
		  fmpz_randbits(rden, state, mn_bits);
		  fmpq_set_fmpz_frac(&mn[k], rnum, rden);
		}
	      
	      flint_printf("delta = %wd; ell = %wd; beta = ", delta, ell);
	      fmpz_poly_print_pretty(beta, "x");
	      flint_printf(", parameters are\n");		      
	      fmpq_print(&mn[0]); flint_printf("\n");
	      fmpq_print(&mn[1]); flint_printf("\n");
	      flint_printf("Gundlach invariants are\n");
	      gundlach_from_hilbert_param(mn, mn, 5);	      
	      fmpq_print(&mn[0]); flint_printf("\n");
	      fmpq_print(&mn[1]); flint_printf("\n");	      
	      res = hilbert_modeq_nonsym_gundlach_eval_Q(num_vec, den, mn, ell, beta, delta);
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

      for (k = 0; k < 2; k++) fmpz_poly_clear(&num_vec[k]);
      fmpz_clear(den);
      fmpz_poly_clear(beta);
      _fmpq_vec_clear(mn, 2);
      fmpz_clear(rnum);
      fmpz_clear(rden);
    }


  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

