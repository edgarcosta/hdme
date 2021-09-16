
#include "igusa.h"

int thomae_keep_candidate(const acb_mat_t tau, acb_srcptr I, slong prec)
{
  acb_ptr j;
  acb_ptr j_test;
  int aux, res;
  slong k;

  j = _acb_vec_init(3);
  j_test = _acb_vec_init(3);

  if (siegel_not_in_fundamental_domain(tau, prec))
    {
      res = 0;
    }
  else
    {
      aux = igusa_from_tau(j, tau, prec);
      igusa_from_cov(j_test, I, prec);
      
      res = 1;
      if (aux) /* We can say nothing if computing j fails */
	{
	  for (k = 0; k < 3; k++)
	    {
	      if (!acb_overlaps(&j[k], &j_test[k])) res = 0;
	    }
	}
    }
  
  _acb_vec_clear(j, 3);
  _acb_vec_clear(j_test, 3);
  return res;
}
