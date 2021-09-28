
#include "hilbert.h"

int hilbert_linear_combination(fmpz* abcde, const acb_mat_t tau, slong delta, slong prec)
{
  fmpz_lll_t fl;
  fmpz_mat_t B;
  fmpz_mat_t U;
  acb_t coeff;
  slong k;
  int res;

  fmpz_lll_context_init_default(fl);
  fmpz_mat_init(B, 7, 5);
  fmpz_mat_init(U, 5, 5);
  acb_init(coeff);

  /* Set up LLL matrix */
  fmpz_one(fmpz_mat_entry(B, 0, 0));
  fmpz_one(fmpz_mat_entry(B, 1, 1));
  fmpz_one(fmpz_mat_entry(B, 2, 2));
  fmpz_one(fmpz_mat_entry(B, 3, 3));
  fmpz_one(fmpz_mat_entry(B, 4, 4));

  acb_mul_2exp_si(coeff, acb_mat_entry(tau, 0, 0), prec/2);
  arf_get_fmpz(fmpz_mat_entry(B, 5, 0),
	       arb_midref(acb_realref(coeff)),
	       ARF_RND_NEAR);
  arf_get_fmpz(fmpz_mat_entry(B, 6, 0),
	       arb_midref(acb_imagref(coeff)),
	       ARF_RND_NEAR);
  
  acb_mul_2exp_si(coeff, acb_mat_entry(tau, 1, 1), prec/2);
  arf_get_fmpz(fmpz_mat_entry(B, 5, 1),
	       arb_midref(acb_realref(coeff)),
	       ARF_RND_NEAR);  
  arf_get_fmpz(fmpz_mat_entry(B, 6, 1),
	       arb_midref(acb_imagref(coeff)),
	       ARF_RND_NEAR);
  
  acb_mul_2exp_si(coeff, acb_mat_entry(tau, 0, 1), prec/2);
  arf_get_fmpz(fmpz_mat_entry(B, 5, 2),
	       arb_midref(acb_realref(coeff)),
	       ARF_RND_NEAR);  
  arf_get_fmpz(fmpz_mat_entry(B, 6, 2),
	       arb_midref(acb_imagref(coeff)),
	       ARF_RND_NEAR);

  acb_mat_det(coeff, tau, prec);
  acb_mul_2exp_si(coeff, coeff, prec/2);
  arf_get_fmpz(fmpz_mat_entry(B, 5, 3),
	       arb_midref(acb_realref(coeff)),
	       ARF_RND_NEAR);  
  arf_get_fmpz(fmpz_mat_entry(B, 6, 3),
	       arb_midref(acb_imagref(coeff)),
	       ARF_RND_NEAR);

  fmpz_one(fmpz_mat_entry(B, 5, 4));
  fmpz_mul_2exp(fmpz_mat_entry(B, 5, 4),
		fmpz_mat_entry(B, 5, 4),
		prec/2);

  /* Call LLL */
  fmpz_lll(B, U, fl);
  
  flint_printf("(hilbert_linear_combination) The following matrix should a small column:\n");
  fmpz_mat_print_pretty(B);
  flint_printf("\n");
  
  /* Get abcde from first column */
  for (k = 0; k < 5; k++)
    {
      fmpz_set(&abcde[k], fmpz_mat_entry(B, k, 0));
    }

  /* Do some kind of check */
  res = 1;
  
  fmpz_mat_clear(B);
  fmpz_mat_clear(U);
  acb_init(coeff);
  return res;
}
