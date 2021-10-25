
#include "hilbert.h"

void cov_from_gundlach(acb_ptr G, acb_srcptr g, slong delta, slong prec)
{
  acb_ptr res;
  res = _acb_vec_init(3);

  acb_one(&res[0]);
  acb_inv(&res[2], &g[0], prec);
  acb_mul(&res[1], &g[1], &G[2], prec);

  _acb_vec_set(G, res, 3);
  _acb_vec_clear(res, 3);
}
