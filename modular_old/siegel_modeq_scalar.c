
#include "modular.h"

void siegel_modeq_scalar(acb_t scal, acb_srcptr I_tau, acb_srcptr stardets,
			 slong ell, slong prec)
{
  slong e, a, b;
  slong wl = 20 * siegel_nb_cosets(ell);
  acb_t temp, res;
  slong k;

  acb_init(temp);
  acb_init(res);
  
  siegel_modeq_exps(&e, &a, &b, ell);

  acb_one(res);
  
  acb_pow_si(temp, &I_tau[1], e, prec);
  acb_mul(res, res, temp, prec);

  acb_pow_si(temp, &I_tau[3], a, prec);
  acb_div(res, res, temp, prec);
  acb_pow_si(temp, &I_tau[1], b, prec);
  acb_div(res, res, temp, prec);

  acb_set_si(temp, 2);
  acb_pow_si(temp, temp, 15 * (wl/12), prec),
  acb_mul(res, res, temp, prec);
  acb_set_si(temp, 3);
  acb_pow_si(temp, temp, wl/4, prec);
  acb_mul(res, res, temp, prec);

  acb_set_si(temp, 2);
  acb_pow_si(temp, temp, -wl, prec);
  acb_mul(res, res, temp, prec);

  for (k = 0; k < siegel_nb_cosets(ell); k++)
    {
      acb_pow_si(temp, &stardets[k], -20, prec);
      acb_mul(res, res, temp, prec);
    }

  acb_set(scal, res);
  
  acb_clear(temp);
  acb_clear(res);
}
