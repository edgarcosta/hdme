
#include "igusa.h"

void cov_from_theta2(acb_ptr I, acb_srcptr theta2, slong prec)
{
  acb_ptr h;
  h = _acb_vec_init(4);
  
  igusa_h(h, theta2, prec);
  acb_div(&I[0], &h[3], &h[2], prec);
  acb_set(&I[1], &h[0]);
  acb_set(&I[2], &h[1]);
  acb_set(&I[3], &h[2]);

  _acb_vec_clear(h, 4);
}
