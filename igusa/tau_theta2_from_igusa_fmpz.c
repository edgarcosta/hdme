
#include "igusa.h"

int tau_theta2_from_igusa_fmpz(acb_mat_t tau, acb_ptr th2, fmpz* I, slong prec)
{
  acb_ptr I_acb;
  slong k;
  int res;

  I_acb = _acb_vec_init(4);
  for (k = 0; k < 4; k++) acb_set_fmpz(&I_acb[k], &I[k]);

  if (cov_is_g2_curve_fmpz(I))
    {
      res = tau_theta2_from_igusa(tau, th2, I_acb, prec);
    }
  else
    {
      res = tau_theta2_from_igusa_ec(tau, th2, I_acb, prec);
    }
  
  _acb_vec_clear(I, 4);  
  return res;  
}
