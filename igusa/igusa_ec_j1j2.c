
#include "igusa.h"

/* Cf. Igusa, "On Siegel modular forms of genus two", p. 181 */

static void igusa_y(acb_t y1, acb_t y2, acb_srcptr I, slong prec)
{
  acb_pow_si(y1, cov_I4(I), 3, prec);
  acb_div(y1, y1, cov_I12(I), prec);
  acb_mul_si(y1, y1, n_pow(2,11) * 3, prec);

  acb_pow_si(y2, cov_I6prime(I), 2, prec);
  acb_div(y2, y2, cov_I12(I), prec);
  acb_mul_si(y2, y2, n_pow(2,14), prec);
}

void igusa_ec_j1j2(acb_ptr j, acb_srcptr I, slong prec)
{
  acb_poly_t pol;
  acb_t y1, y2;

  acb_poly_init(pol);
  acb_init(y1);
  acb_init(y2);

  igusa_y(y1, y2, I, prec);

  acb_poly_set_coeff_si(pol, 2, 1);
  acb_poly_set_coeff_acb(pol, 0, y1);

  acb_sub(y2, y2, y1, prec);
  acb_sub_si(y2, y2, n_pow(2,12) * n_pow(3,6), prec);
  acb_div_si(y2, y2, n_pow(2,6) * n_pow(3,3), prec);
  acb_poly_set_coeff_acb(pol, 1, y2);

  /* Assume that thomae_roots will not throw for deg 2 polynomials */
  thomae_roots(j, pol, prec);

  acb_poly_clear(pol);
  acb_clear(y1);
  acb_clear(y2);
}
