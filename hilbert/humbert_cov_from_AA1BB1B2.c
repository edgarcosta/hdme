
#include "hilbert.h"

void humbert_cov_from_AA1BB1B2(acb_ptr I, acb_srcptr AA1BB1B2, slong prec)
{
  acb_t temp;
  acb_ptr res;

  acb_init(temp);
  res = _acb_vec_init(4);

  acb_mul_si(temp, &AA1BB1B2[3], -24, prec);
  acb_div(&res[0], temp, &AA1BB1B2[1], prec);

  acb_mul_si(&res[1], &AA1BB1B2[0], -12, prec);

  /* Cf formula for I6prime */
  acb_mul_si(&res[2], &AA1BB1B2[2], 54, prec);

  acb_mul_si(temp, &AA1BB1B2[1], -4, prec);
  acb_mul(&res[2], temp, &AA1BB1B2[4], prec);

  _acb_vec_set(I, res, 4);
  acb_clear(temp);
  _acb_vec_clear(res, 4);
}
