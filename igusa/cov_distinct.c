
#include "igusa.h"

int cov_distinct(acb_srcptr I1, acb_srcptr I2, slong prec)
{
  acb_ptr test;
  slong j;
  int r1, r2;

  test = _acb_vec_init(4);

  for (j = 0; j < 4; j++) acb_div(&test[j], &I1[j], &I2[j], prec);
  r1 = cov_no_rescale_to_one(test);  
  for (j = 0; j < 4; j++) acb_div(&test[j], &I2[j], &I1[j], prec);
  r2 = cov_no_rescale_to_one(test);

  _acb_vec_clear(test, 4);
  return (r1 || r2);
}
