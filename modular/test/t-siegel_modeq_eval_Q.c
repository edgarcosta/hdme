
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
      fmpz_poly_struct num_vec[3];
      fmpz_t den;
      fmpq* j;
      slong ell =  2; /* n_randprime(state, 3, 1); */ /* 2, 3, 5 or 7 */
      slong j_bits = 1 + n_randint(state, 20);
      slong k;
      int res;

      for (k = 0; k < 3; k++) fmpz_poly_init(&num_vec[k]);
      fmpz_init(den);
      j = _fmpq_vec_init(3);

      for (k = 0; k < 3; k++)
	{
	  fmpq_randtest_not_zero(&j[k], state, j_bits);
	}
      res = siegel_modeq_eval_Q(num_vec, den, j, ell);
      if (!res)
	{
	  flint_printf("FAIL\n");
	  for (k = 0; k < 3; k++)
	    {
	      fmpq_print(&j[k]); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}
      
      for (k = 0; k < 3; k++) fmpz_poly_clear(&num_vec[k]);
      fmpz_clear(den);
      _fmpq_vec_clear(j, 3);
    }
  
  
  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
