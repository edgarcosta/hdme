
#include "igusa.h"

int igusa_from_tau(acb_ptr j, const acb_mat_t tau, slong prec)
{
  slong g = 2;
  acb_ptr th2;
  int res;
  
  th2 = _acb_vec_init(n_pow(2, 2*g));
  
  res = theta2_unif(th2, tau, prec);
  igusa_from_theta2(j, th2, prec);
  
  _acb_vec_clear(th2, n_pow(2, 2*g));
  return res;
}
