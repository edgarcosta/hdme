
#include "igusa.h"

void cov_from_h(acb_ptr I, acb_srcptr h, slong prec)
{
  acb_ptr res;
  res = _acb_vec_init(4);
  
  acb_div(&res[0], &h[3], &h[2], prec);
  acb_set(&res[1], &h[0]);
  acb_set(&res[2], &h[1]);
  acb_set(&res[3], &h[2]);

  _acb_vec_set(I, res, 4);
  _acb_vec_clear(res, 4);
}
