
#include "igusa.h"

void cov_from_theta2(acb_ptr I, acb_srcptr theta2, slong prec)
{
  /* Support weird aliasing... */
  acb_ptr res;
  res = _acb_vec_init(4);
  
  igusa_h4(cov_I4(res), theta2, prec);
  igusa_h6(cov_I6prime(res), theta2, prec);
  igusa_h10(cov_I10(res), theta2, prec);
  igusa_h12(cov_I12(res), theta2, prec);

  _acb_vec_set(I, res, 4);
  _acb_vec_clear(res, 4);
}
