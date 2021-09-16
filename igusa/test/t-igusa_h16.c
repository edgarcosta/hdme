
#include "igusa.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("igusa_h....");
  fflush(stdout);

  flint_randinit(state);

  /* Test formula I6 = h16/h10 */
  for (iter = 0; iter < 1000 * arb_test_multiplier(); iter++)
    {
      acb_ptr theta2;
      acb_ptr h, I;
      acb_t h16;
      acb_t I6;
      slong k;

      slong prec = 50 + n_randint(state, 2000);

      theta2 = _acb_vec_init(16);
      h = _acb_vec_init(4);
      I = _acb_vec_init(4);
      acb_init(h16);
      acb_init(I6);

      theta2_randtest(theta2, state, prec);
      igusa_h(h, theta2, prec);
      
      acb_div(&I[0], &h[3], &h[2], prec);
      acb_set(&I[1], &h[0]);
      acb_set(&I[2], &h[1]);
      acb_set(&I[3], &h[2]);

      igusa_I6(I6, I, prec);
      igusa_h16(h16, theta2, prec);
      acb_div(h16, h16, &h[2], prec);

      if (!acb_overlaps(h16, I6))
	{
	  flint_printf("FAIL (h16)");
	  flint_printf("Theta constants:\n");
	  for (k = 0; k < 16; k++)
	    {
	      acb_printd(&theta2[k], 30); flint_printf("\n");
	    }
	  flint_printf("\nh4: "); acb_printd(&h[0], 30);
	  flint_printf("\nh6: "); acb_printd(&h[1], 30);
	  flint_printf("\nh10: "); acb_printd(&h[2], 30);
	  flint_printf("\nh12: "); acb_printd(&h[3], 30);
	  flint_printf("\nh16/h10: "); acb_printd(h16, 30);
	  flint_printf("\nI6: "); acb_printd(I6, 30);
	  fflush(stdout);
	  flint_abort();
	}

      _acb_vec_clear(theta2, 16);
      _acb_vec_clear(h, 4);
      _acb_vec_clear(I, 4);
      acb_clear(h16);
      acb_clear(I6);
    }
 
  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

      

