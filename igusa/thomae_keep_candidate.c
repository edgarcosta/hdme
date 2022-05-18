
#include "igusa.h"

int thomae_keep_candidate(const acb_mat_t tau, acb_srcptr I, slong prec)
{
  acb_ptr test;
  int aux, res;

  test = _acb_vec_init(4);

  if (siegel_not_in_fundamental_domain(tau, prec))
    {
      res = 0;
    }
  else
    {
      aux = cov_from_tau(test, tau, prec);      
      res = 1;
      if (aux && cov_distinct(test, I, prec)) res = 0;
      /* We can say nothing if computing test fails */
    }
  
  _acb_vec_clear(test, 4);
  return res;
}
