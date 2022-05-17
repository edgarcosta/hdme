
#include "igusa.h"

void cov_rescale(acb_ptr I, acb_srcptr S, const acb_t scal, slong prec)
{
  acb_t temp;
  acb_ptr res;
  slong wt[4] = COV_WEIGHTS;
  slong j;
  
  acb_init(temp);
  res = _acb_vec_init(4);

  for (j = 0; j < 4; j++)
    {
      acb_pow_ui(temp, scal, wt[j], prec);
      acb_mul(&res[j], &S[j], temp, prec);
    }
  
  _acb_vec_set(I, res, 4);
  _acb_vec_clear(res, 4);
  acb_clear(temp);
}
