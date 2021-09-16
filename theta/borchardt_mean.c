
#include "theta.h"

int borchardt_mean(acb_t r, acb_srcptr b, slong prec)
{
  acb_ptr a;
  fmpz_t nb;
  slong i;
  int res;

  a = _acb_vec_init(4);
  fmpz_init(nb);

  /* Set a to (1, b1/b0, b2/b0, b3/b0) */
  _acb_vec_scalar_div(a, b, 4, &b[0], prec);
  acb_one(&a[0]);

  /* Adjust target precision given input */
  for (i = 0; i < 4; i++)
    {
      prec = FLINT_MIN(prec, -acb_rel_error_bits(&a[i]));
      prec = FLINT_MAX(prec, 10);
    }
  
  res = borchardt_mean_nb_steps_before_quad_conv(nb, a, prec);

  /* flint_printf("\nprec = %wd\n", prec);
     flint_printf("Steps before quadratic convergence: "); fmpz_print(nb); flint_printf("\n"); */
  i = 0;

  if (res)
    {
      for (i = 0; fmpz_cmp_si(nb, i) > 0; i++)
	{
	  res = borchardt_step(a, a, prec);
	  if (!res) break;
	  if (borchardt_mean_quad_conv_is_reached(a, prec)) break;
	}
    }
  
  /* flint_printf("Quadratic convergence after %wd steps\n", i+1); */
  if (res) borchardt_mean_nb_steps_after_quad_conv(nb, a, prec);
  /* for (i = 0; i < 4; i++)
     {
     flint_printf("a[%wd] = ", i); acb_printd(&a[i], prec); flint_printf("\n");
     }
     flint_printf("Steps after quadratic convergence: "); fmpz_print(nb); flint_printf("\n"); */
  
  /* After nb Borchardt steps, the remaining error is at most 2^(-prec) */
  if (res)
    {
      for (i = 0; fmpz_cmp_si(nb, i) > 0; i++)
	{
	  res = borchardt_step(a, a, prec);
	  if (!res) break;
	}

      /* flint_printf("\nBorchardt result:\n");
	 for (i = 0; i < 4; i++)
	 {
	 flint_printf("a[%wd] = ", i); acb_printd(&a[i], prec); flint_printf("\n");
	 }
      */
      
      acb_set(r, &a[0]);
      arb_add_error_2exp_si(acb_realref(r), -prec);
      arb_add_error_2exp_si(acb_imagref(r), -prec);
      acb_mul(r, r, &b[0], prec);
    }
  
  _acb_vec_clear(a, 4);
  fmpz_clear(nb);
  return res;
}
