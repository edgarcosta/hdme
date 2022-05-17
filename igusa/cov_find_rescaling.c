
#include "igusa.h"

int cov_find_rescaling(acb_t scal, acb_srcptr I, fmpz* S, slong prec)
{
  slong wt, i1, i2, e1, e2;
  acb_t r1, r2;
  acb_t ci, cs;
  acb_ptr rescale;
  slong j;
  int res = 1;

  acb_init(r1);
  acb_init(r2);
  acb_init(ci);
  acb_init(cs);
  rescale = _acb_vec_init(4);

  /* Compute combination of S with minimum weight */
  cov_min_weight_combination(&wt, &i1, &i2, &e1, &e2, S);
  acb_set_fmpz(r1, &S[i1]);
  acb_pow_si(r1, r1, e1, prec);
  acb_set_fmpz(r2, &S[i2]);
  acb_pow_si(r2, r2, e2, prec);
  acb_mul(cs, r1, r2, prec);

  /* Repeat same combination for I, quotient and root */
  acb_pow_si(r1, &I[i1], e1, prec);
  acb_pow_si(r2, &I[i2], e2, prec);
  acb_mul(ci, r1, r2, prec);
  acb_div(cs, cs, ci, prec);
  borchardt_root_ui(cs, cs, wt, prec);

  /* Check that rescaling works */
  cov_rescale(rescale, I, cs, prec);
  for (j = 0; j < 4; j++)
    {
      if (!acb_contains_fmpz(&rescale[j], &S[j]))
	{
	  res = 0;
	  break;
	}
    }
  
  return res;
}
