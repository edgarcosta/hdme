
#include "modular.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("modeq_rational_poly....");
  fflush(stdout);
  
  flint_randinit(state);
  
  for (iter = 0; iter < 100 * arb_test_multiplier(); iter++)
    {
      fmpq_t c;
      fmpq_t r;
      fmpz_t den;
      acb_t x;
      slong prec = 100 + n_randint(state, n_pow(10, 4));
      slong bits = prec/4;
      int res;

      fmpq_init(c);
      fmpq_init(r);
      fmpz_init(den);
      acb_init(x);
      
      fmpq_randtest(c, state, bits);
      fmpz_one(den);
      acb_set_fmpq(x, c, prec);
      res = modeq_rational_coeff(r, den, x, den, prec);
      if (!res || !fmpq_equal(r, c))
	{
	  flint_printf("FAIL (prec = %wd)\n", prec);
	  fmpq_print(c); flint_printf("\n");
	  acb_printd(x, 30); flint_printf("\n");
	  fmpq_print(r); flint_printf("\n");	  
	  fflush(stdout);
	  flint_abort();
	}

      fmpq_randbits(c, state, 2 * prec/3);
      acb_set_fmpq(x, c, prec);
      res = modeq_rational_coeff(r, den, x, den, prec);
      if (res && !fmpq_equal(r, c))
	{
	  flint_printf("FAIL (spurious coefficient recognized)\n", prec);
	  fmpq_print(c); flint_printf("\n");
	  acb_printd(x, 30); flint_printf("\n");
	  fmpq_print(r); flint_printf("\n");	  
	  fflush(stdout);
	  flint_abort();
	}
      
      fmpq_clear(c);
      fmpq_clear(r);
      fmpz_clear(den);
      acb_clear(x);
    }
  
  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

