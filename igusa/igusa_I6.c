
#include "igusa.h"

void igusa_I6(acb_t I6, acb_srcptr I, slong prec)
{
  acb_t I2, res, temp;

  acb_init(I2);
  acb_init(res);
  acb_init(temp);

  acb_div(I2, igusa_I12(I), igusa_I10(I), prec);
  acb_mul_si(res, igusa_I6prime(I), 2, prec);
  acb_mul(temp, I2, igusa_I4(I), prec);
  acb_sub(res, res, temp, prec);
  acb_div_si(res, res, -3, prec);

  acb_set(I6, res);

  acb_clear(I2);
  acb_clear(res);
  acb_clear(temp);
}
