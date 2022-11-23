
#include "theta.h"

int main()
{
  slong iter;
  flint_rand_t state;

  flint_printf("borchardt_mean_invalid....");
  fflush(stdout);

  flint_randinit(state);

  for (iter = 0; iter < 500 * arb_test_multiplier(); iter++)
    {
      acb_ptr a;
      slong prec = 50 + n_randint(state, 500);
      slong mag_bits = 1 + n_randint(state, 5);
      int res;
      slong k;

      a = _acb_vec_init(4);

      /* If all entries have positive or negative real part, then it's
         not invalid */
      for (k = 0; k < 4; k++)
        {
          acb_randtest_precise(&a[k], state, prec, mag_bits);
          arb_sqr(acb_realref(&a[k]), acb_realref(&a[k]), prec);
        }
      res = borchardt_mean_invalid(a, prec);
      if (res)
        {
          flint_printf("FAIL (positive real parts)\n");
          for (k = 0; k < 4; k++)
            {
              acb_printd(&a[k], 30); flint_printf("\n");
            }
          fflush(stdout);
          flint_abort();
        }

      for (k = 0; k < 4; k++)
        {
          acb_randtest_precise(&a[k], state, prec, mag_bits);
          arb_sqr(acb_realref(&a[k]), acb_realref(&a[k]), prec);
          acb_neg(&a[k], &a[k]);
        }
      res = borchardt_mean_invalid(a, prec);
      if (res)
        {
          flint_printf("FAIL (negative real parts)\n");
          for (k = 0; k < 4; k++)
            {
              acb_printd(&a[k], 30); flint_printf("\n");
            }
          fflush(stdout);
          flint_abort();
        }


      /* If there is one entry on each half axis, then it is invalid */
      for (k = 0; k < 4; k++) acb_zero(&a[k]);
      arf_randtest_not_zero(arb_midref(acb_realref(&a[0])), state, prec, mag_bits);
      arb_sqr(acb_realref(&a[0]), acb_realref(&a[0]), prec);

      arf_randtest_not_zero(arb_midref(acb_imagref(&a[1])), state, prec, mag_bits);
      arb_sqr(acb_imagref(&a[1]), acb_imagref(&a[1]), prec);

      arf_randtest_not_zero(arb_midref(acb_realref(&a[2])), state, prec, mag_bits);
      arb_sqr(acb_realref(&a[2]), acb_realref(&a[2]), prec);
      acb_neg(&a[2], &a[2]);

      arf_randtest_not_zero(arb_midref(acb_imagref(&a[3])), state, prec, mag_bits);
      arb_sqr(acb_imagref(&a[3]), acb_imagref(&a[3]), prec);
      acb_neg(&a[3], &a[3]);

      res = borchardt_mean_invalid(a, prec);
      if (!res)
        {
          flint_printf("FAIL (axes)\n");
          for (k = 0; k < 4; k++)
            {
              acb_printd(&a[k], 30); flint_printf("\n");
            }
          fflush(stdout);
          flint_abort();
        }

      _acb_vec_clear(a, 4);
    }
  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
