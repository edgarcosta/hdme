
#include "hecke.h"

int hecke_set_I(hecke_t H, fmpz* I, slong prec)
{
  int res;
  slong k;
  
  for (k = 0; k < 4; k++)
    {
      acb_set_fmpz(&hecke_I_tau(H)[k], &I[k]);
    }
  res = tau_theta2_from_igusa(hecke_tau(H), hecke_theta2_tau(H),
			      hecke_I_tau(H), prec);
  return res;
}
