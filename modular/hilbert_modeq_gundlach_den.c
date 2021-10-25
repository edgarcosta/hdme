
#include "modular.h"

void hilbert_modeq_gundach_den(acb_t den, acb_srcptr I_vec_beta,
			       acb_srcptr I_vec_betabar, const acb_t scal,
			       slong ell, slong delta, slong prec)
{
  acb_t res;
  acb_ptr G;
  slong k;

  acb_init(res);
  G = _acb_vec_init(3);
  
  acb_one(res);
  for (k = 0; k < hilbert_nb_cosets(ell, delta); k++)
    {
      gundlach_cov_from_igusa(G, &I_vec_beta[4*k], delta, prec);
      acb_mul(res, res, &G[2], prec);
    }
  acb_mul(res, res, scal, prec);
  acb_set(den, res);

  acb_clear(res);
  _acb_vec_clear(G, 3);
}
