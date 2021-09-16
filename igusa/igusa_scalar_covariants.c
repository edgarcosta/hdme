
#include "igusa.h"

void igusa_scalar_covariants(acb_ptr I, const acb_poly_t crv, slong prec)
{
  /* Support weird aliasing */
  acb_ptr res;
  acb_t a0, a1, a2, a3, a4, a5, a6;
  res = _acb_vec_init(4);

  acb_init(a0);
  acb_init(a1);
  acb_init(a2);
  acb_init(a3);
  acb_init(a4);
  acb_init(a5);
  acb_init(a6);

  acb_poly_get_coeff_acb(a0, crv, 0);
  acb_poly_get_coeff_acb(a1, crv, 1);
  acb_poly_get_coeff_acb(a2, crv, 2);
  acb_poly_get_coeff_acb(a3, crv, 3);
  acb_poly_get_coeff_acb(a4, crv, 4);
  acb_poly_get_coeff_acb(a5, crv, 5);
  acb_poly_get_coeff_acb(a6, crv, 6);

  igusa_I2_autogen(&res[0], a0, a1, a2, a3, a4, a5, a6, prec);
  igusa_I4_autogen(&res[1], a0, a1, a2, a3, a4, a5, a6, prec);
  igusa_I6prime_autogen(&res[2], a0, a1, a2, a3, a4, a5, a6, prec);
  igusa_I10_autogen(&res[3], a0, a1, a2, a3, a4, a5, a6, prec);
  _acb_vec_set(I, res, 4);
  
  _acb_vec_clear(res, 4);
  acb_clear(a0);
  acb_clear(a1);
  acb_clear(a2);
  acb_clear(a3);
  acb_clear(a4);
  acb_clear(a5);
  acb_clear(a6);
}
