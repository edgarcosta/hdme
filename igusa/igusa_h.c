
#include "igusa.h"

void igusa_h(acb_ptr h, acb_srcptr theta2, slong prec)
{
  /* Support weird aliasing... */
  acb_ptr res;
  res = _acb_vec_init(4);
  
  igusa_h4(&res[0], theta2, prec);
  igusa_h6(&res[1], theta2, prec);
  igusa_h10(&res[2], theta2, prec);
  igusa_h12(&res[3], theta2, prec);

  _acb_vec_set(h, res, 4);
  _acb_vec_clear(res, 4);
}
