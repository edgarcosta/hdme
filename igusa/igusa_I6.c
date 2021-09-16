
#include "igusa.h"

void igusa_I6(acb_t I6, acb_srcptr I, slong prec)
{
  acb_t res, temp;

  acb_init(res);
  acb_init(temp);
  
  /* Get I6 from I6prime */
  acb_mul_si(res, &I[2], 2, prec);
  acb_mul(temp, &I[0], &I[1], prec);
  acb_sub(res, res, temp, prec);
  acb_div_si(res, res, -3, prec);

  acb_set(I6, res);
  acb_clear(res);
  acb_clear(temp);
}
