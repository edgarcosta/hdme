
#include "igusa.h"

int main()
{
  
  slong iter;
  flint_rand_t state;
  
  flint_printf("mestre....");
  fflush(stdout);

  flint_randinit(state);

  for (iter = 0; iter < 50 * arb_test_multiplier(); iter++)
    {
      slong prec = 500 + n_randint(state, 1000);
      acb_poly_t crv_test, crv;
      acb_ptr I_test, I;
      acb_ptr j_test, j;
      int res;
      int k;

      acb_poly_init(crv_test);
      acb_poly_init(crv);
      I_test =  _acb_vec_init(4);
      I =  _acb_vec_init(4);
      j_test =  _acb_vec_init(3);
      j =  _acb_vec_init(3);

      igusa_generic_randtest(crv_test, I_test, state, prec);
      igusa_from_curve(j_test, crv_test, prec);
      res = igusa_is_defined(j_test);

      if (!res)
	{
	  flint_printf("FAIL (j_test)\n");
	  flint_printf("Curve: "); acb_poly_printd(crv_test, 30); flint_printf("\n");
	  fflush(stdout);
	  flint_abort();
	}
      
      mestre(crv, I_test, prec);
      igusa_from_curve(j, crv, prec);
      res = igusa_is_defined(j);
      
      if (!res)
	{
	  flint_printf("FAIL (j)\n");
	  flint_printf("Curve: "); acb_poly_printd(crv, 30); flint_printf("\n");
	  fflush(stdout);
	  flint_abort();
	}
      
      if (!acb_overlaps(&j[0], &j_test[0])
	  || !acb_overlaps(&j[1], &j_test[1])
	  || !acb_overlaps(&j[2], &j_test[2]) )
	{
	  flint_printf("FAIL (overlap)\n");
	  flint_printf("crv_test: "); acb_poly_printd(crv_test, 30); flint_printf("\n");
	  flint_printf("crv: "); acb_poly_printd(crv, 30); flint_printf("\n");
	  for (k = 0; k < 3; k++)
	    {
	      flint_printf("j_test[%d]: ", k); acb_printd(&j_test[k], 30); flint_printf("\n");
	      flint_printf("j[%d]: ", k); acb_printd(&j[k], 30); flint_printf("\n");
	    }
	  fflush(stdout);
	  flint_abort();
	}

      acb_poly_clear(crv_test);
      acb_poly_clear(crv);
      _acb_vec_clear(I_test, 4);
      _acb_vec_clear(I, 4);
      _acb_vec_clear(j_test, 3);
      _acb_vec_clear(j, 3);
    }
  
  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
