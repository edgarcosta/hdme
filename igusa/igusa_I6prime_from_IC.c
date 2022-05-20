
#include "igusa.h"

void igusa_I6prime_from_IC(acb_t I6prime, acb_srcptr I, slong prec)
{
  acb_t res, temp;

  acb_init(res);
  acb_init(temp);
  
  /* Get I6prime from I6 */
  acb_mul_si(res, &I[2], -3, prec);
  acb_mul(temp, &I[0], &I[1], prec);
  acb_add(res, res, temp, prec);
  acb_div_si(res, res, 2, prec);

  acb_set(I6prime, res);
  acb_clear(res);
  acb_clear(temp);
}
