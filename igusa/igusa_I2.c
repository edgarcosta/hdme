
#include "igusa.h"

void igusa_I2(acb_t I2, acb_srcptr I, slong prec)
{
  acb_div(I2, cov_I12(I), cov_I10(I), prec),
}
