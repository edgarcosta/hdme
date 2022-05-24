
#include "modular.h"

void siegel_modeq_scalar(acb_t c, const hecke_t H, fmpz* I,
			 const modeq_ctx_t ctx, slong prec)
{
  slong wt = modeq_ctx_weight(ctx);
  slong ell = hecke_ell(H);
  slong weights[4] = IGUSA_WEIGHTS;
  acb_t s;
  slong k;

  acb_init(s);
  
  /* Step 1: product of all stardets to the correct weight, to get
     a Siegel modular form */
  acb_one(c);
  for (k = 0; k < hecke_nb(H); k++)
    {
      acb_inv(s, hecke_stardet(H, k), prec);
      acb_mul_si(s, s, ell*ell, prec);
      acb_mul(c, c, s, prec);
    }
  
  /* Step 2: adjust power of ell that keeps integer coefficients */
  acb_set_si(s, ell);
  acb_pow_si(s, s, - (2*ell*ell + ell + 1), prec);
  acb_mul(c, c, s, prec);
  
  /* Step 3: multiply by scaling factor between hecke_I_tau and I */
  cov_find_rescaling(s, hecke_I_tau(H), I, 4, weights, prec);  
  acb_mul(c, c, s, prec);

  /* Step 4: compute power according to weight */
  acb_pow_si(c, c, wt, prec);

  /* Step 5: additional factor according to additional MF generators? */

  acb_clear(s);
}
