
#include "igusa.h"

int main()
{
  
  slong iter;
  flint_rand_t state;
  
  flint_printf("cov_from_igusa_fmpz....");
  fflush(stdout);

  flint_randinit(state);

  for (iter = 0; iter < 500 * arb_test_multiplier(); iter++)
    {
      fmpq* j_test;
      fmpz* I;
      fmpq* j;
      fmpq_t scal;
      slong k;
      slong bits = 10 + n_randint(state, 50);
      int res = 1;

      j_test = _fmpq_vec_init(3);
      I = _fmpz_vec_init(4);
      j = _fmpq_vec_init(3);
      fmpq_init(scal);

      fmpq_randtest_not_zero(scal, state, bits);
      for (k = 0; k < 3; k++)
	{
	  fmpq_randtest_not_zero(&j[k], state, bits);
	  fmpq_mul(&j[k], &j[k], scal);
	}
      cov_from_igusa_fmpz(I, j);
      igusa_from_cov_fmpz(j_test, I);

      for (k = 0; k < 3; k++)
	{
	  if (!fmpq_equal(&j[k], &j_test[k])) res = 0;
	}
      if (!res)
	{
	  flint_printf("FAIL (j != j_test)\n");
	  for (k = 0; k < 3; k++)
	    {
	      fmpq_print(&j_test[k]);
	      flint_printf("\n");
	      fmpq_print(&j[k]);
	      flint_printf("\n");
	      fflush(stdout);
	    }
	  flint_abort();
	}

      _fmpq_vec_clear(j_test, 3);
      _fmpz_vec_clear(I, 4);
      _fmpq_vec_clear(j, 3);
      fmpq_clear(scal);
    }
  
  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

