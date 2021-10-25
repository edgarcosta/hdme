
#include "igusa.h"

void cov_rescale(acb_ptr I, acb_srcptr S, const acb_t scal, slong prec)
{
  acb_t temp;
  acb_ptr res;
  
  acb_init(temp);
  res = _acb_vec_init(4);

  acb_pow_ui(temp, scal, 2, prec);
  acb_mul(&res[0], &S[0], temp, prec);
  
  acb_pow_ui(temp, scal, 4, prec);
  acb_mul(&res[1], &S[1], temp, prec);
  
  acb_pow_ui(temp, scal, 6, prec);
  acb_mul(&res[2], &S[2], temp, prec);
  
  acb_pow_ui(temp, scal, 10, prec);
  acb_mul(&res[3], &S[3], temp, prec);

  _acb_vec_set(I, res, 4);
  _acb_vec_clear(res, 4);
  acb_clear(temp);
}
