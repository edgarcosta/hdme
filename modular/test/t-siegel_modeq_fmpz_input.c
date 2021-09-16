
#include "modular.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("siegel_modeq_fmpz_input....");
  fflush(stdout);

  flint_randinit(state);

  for (iter = 0; iter < 100 * arb_test_multiplier(); iter++)
    {
      fmpz_t den;
      fmpz* num;
      fmpq* j;
      fmpq* j_test;
      slong len = 1 + n_randint(state, 10);
      slong k;
      slong bits = 1 + n_randint(state, 10);
      int res = 1;

      fmpz_init(den);
      num = _fmpz_vec_init(len);
      j = _fmpq_vec_init(len);
      j_test = _fmpq_vec_init(len);

      for (k = 0; k < len; k++)
	{
	  fmpq_randtest(&j[k], state, bits);
	}
      siegel_modeq_fmpz_input(den, num, j, len);
      for (k = 0; k < len; k++)
	{
	  fmpq_set_fmpz_frac(&j_test[k], &num[k], den);
	  if (!fmpq_equal(&j[k], &j_test[k])) res = 0;
	}
      if (!res)
	{
	  flint_printf("FAIL\n");
	  for (k = 0; k < len; k++)
	    {
	      fmpq_print(&j[k]); flint_printf("\n");
	      fmpq_print(&j_test[k]); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}

      fmpz_clear(den);
      _fmpz_vec_clear(num, len);
      _fmpq_vec_clear(j, len);
      _fmpq_vec_clear(j_test, len);
    }
  
  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
