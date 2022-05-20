
#include "modular.h"

void modeq_cov(acb_ptr I_vec, acb_srcptr th2_vec, slong nb, slong prec)
{
  slong k;
  for (k = 0; k < nb; k++)
    {
      cov_from_theta2(&I_vec[4*k], &th2_vec[16*k], prec);
    }
}
