
#include "igusa.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("igusa_scalar_covariants....");
  fflush(stdout);

  flint_randinit(state);

  for (iter = 0; iter < 100 * arb_test_multiplier(); iter++)
    {
      slong prec = 50 + n_randint(state, 1000);
      acb_ptr roots;
      acb_ptr I;
      acb_t R2;
      acb_poly_t crv;

      /* Generate crv with multiple roots; check that covariants are zero */
      roots = _acb_vec_init(6);
      I = _acb_vec_init(4);
      acb_poly_init(crv);
      acb_init(R2);

      acb_randtest_precise(&roots[0], state, prec, 1);
      acb_set(&roots[1], &roots[0]);
      acb_set(&roots[2], &roots[0]);
      acb_set(&roots[3], &roots[0]);
      acb_randtest_precise(&roots[4], state, prec, 1);
      acb_randtest_precise(&roots[5], state, prec, 1);

      acb_poly_product_roots(crv, roots, 6, prec);

      igusa_scalar_covariants(I, crv, prec);

      if (!acb_contains_zero(&I[0])
	  || !acb_contains_zero(&I[1])
	  || !acb_contains_zero(&I[2])
	  || !acb_contains_zero(&I[3]))
	{
	  flint_printf("FAIL\n");
	  flint_printf("Curve: "), acb_poly_printd(crv, 30), flint_printf("\n");
	  flint_printf("I2 = "); acb_printd(&I[0], 30); flint_printf("\n");
	  flint_printf("I4 = "); acb_printd(&I[1], 30); flint_printf("\n");
	  flint_printf("I6prime = "); acb_printd(&I[2], 30); flint_printf("\n");
	  flint_printf("I10 = "); acb_printd(&I[3], 30); flint_printf("\n");
	  fflush(stdout);
	  flint_abort();
	}
      
      /* Generate crv with less multiple roots; check that R2 is zero */
      acb_randtest_precise(&roots[0], state, prec, 1);
      acb_set(&roots[1], &roots[0]);
      acb_randtest_precise(&roots[2], state, prec, 1);
      acb_set(&roots[3], &roots[2]);
      acb_randtest_precise(&roots[4], state, prec, 1);
      acb_randtest_precise(&roots[5], state, prec, 1);

      acb_poly_product_roots(crv, roots, 6, prec);

      igusa_scalar_covariants(I, crv, prec);
      igusa_I6(R2, I, prec);
      igusa_R2(R2, I, prec);
      
      if (!acb_contains_zero(&I[3])
	  || !acb_contains_zero(R2))
	{
	  flint_printf("FAIL\n");
	  flint_printf("Curve: "), acb_poly_printd(crv, 30), flint_printf("\n");
	  flint_printf("I2 = "); acb_printd(&I[0], 30); flint_printf("\n");
	  flint_printf("I4 = "); acb_printd(&I[1], 30); flint_printf("\n");
	  flint_printf("I6prime = "); acb_printd(&I[2], 30); flint_printf("\n");
	  flint_printf("I10 = "); acb_printd(&I[3], 30); flint_printf("\n");
	  flint_printf("R2 = "); acb_printd(R2, 30); flint_printf("\n");
	  fflush(stdout);
	  flint_abort();
	}
      
      _acb_vec_clear(roots, 6);
      _acb_vec_clear(I, 4);
      acb_poly_clear(crv);
      acb_clear(R2);
    }
  
  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}


      
