
#include "igusa.h"

void cov_from_igusa(acb_ptr I, acb_srcptr j, slong prec)
{
  acb_ptr res;

  res = _acb_vec_init(4);

  acb_set_si(&res[3], 1);
  borchardt_root_ui(&res[1], &j[2], 5, prec);

  acb_div(&res[2], &j[0], &res[1], prec);

  acb_pow_si(&res[0], &res[1], -2, prec);
  acb_mul(&res[0], &res[0], &j[1], prec);

  _acb_vec_set(I, res, 4);
  _acb_vec_clear(res, 4);
}
