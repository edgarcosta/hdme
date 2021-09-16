
#include "igusa.h"

void igusa_from_cov(acb_ptr j, acb_srcptr I, slong prec)
{
  acb_ptr res;
  acb_t I10inv;

  res = _acb_vec_init(3);
  acb_init(I10inv);

  acb_inv(I10inv, &I[3], prec);
  acb_mul(&res[0], &I[1], &I[2], prec);
  acb_mul(&res[0], &res[0], I10inv, prec);
  
  acb_sqr(&res[1], &I[1], prec);
  acb_mul(&res[1], &res[1], &I[0], prec);
  acb_mul(&res[1], &res[1], I10inv, prec);
  
  acb_sqr(I10inv, I10inv, prec); /* Now 1/I10^2 */
  acb_pow_ui(&res[2], &I[1], 5, prec);
  acb_mul(&res[2], &res[2], I10inv, prec);
  
  _acb_vec_set(j, res, 3);
  _acb_vec_clear(res, 3);
  acb_clear(I10inv);
}
