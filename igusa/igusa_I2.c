
#include "igusa.h"

void igusa_I2(acb_t I2, acb_srcptr I, slong prec)
{
  acb_div(I2, igusa_I12(I), igusa_I10(I), prec);
}
