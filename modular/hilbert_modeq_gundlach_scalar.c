
#include "modular.h"

void hilbert_modeq_gundlach_scalar(acb_t scal, acb_srcptr I_tau, acb_srcptr stardets,
				   slong ell, slong delta, slong prec)
{
  slong e, a, b;
  acb_ptr G;
  slong wb = 2 * 10 * hilbert_nb_cosets(ell);
  acb_t temp, res;
  slong k;
  
  if (delta != 5)
    {
      flint_printf("(hilbert_modeq_gundlach_scalar) Error: Gundlach invariants only implemented for discriminant 5\n");
      fflush(stdout);
      flint_abort();
    }

  acb_init(temp);
  acb_init(res);
  G = _acb_vec_init(3);
  
  hilbert_modeq_gundlach_exps(&e, &a, &b, ell);
  gundlach_cov_from_igusa(G, I_tau, delta, prec);

  acb_one(res);
  
  acb_pow_si(temp, &G[0], e, prec);
  acb_mul(res, res, temp, prec);

  acb_pow_si(temp, &G[2], a, prec);
  acb_div(res, res, temp, prec);
  acb_pow_si(temp, &G[1], b, prec);
  acb_div(res, res, temp, prec);

  acb_set_si(temp, 2);
  acb_pow_si(temp, temp, wb, prec);
  acb_mul(res, res, temp, prec); /* We can do better here */

  for (k = 0; k < 2*hilbert_nb_cosets(ell); k++)
    {
      acb_pow_si(temp, &stardets[k], -10, prec);
      acb_mul(res, res, temp, prec);
    }

  acb_set(scal, res);
  
  acb_clear(temp);
  acb_clear(res);
  _acb_vec_clear(G, 3);
}
