
#include "igusa.h"

void igusa_IC(acb_ptr IC, acb_srcptr I, slong prec)
{
  acb_t I2, I6;

  acb_init(I2);
  acb_init(I6);

  igusa_I2(I2, I, prec);
  igusa_I6(I6, I, prec);

  acb_set(&IC[0], I2);
  acb_set(&IC[1], igusa_I4(I));
  acb_set(&IC[2], I6);
  acb_set(&IC[3], igusa_I10(I));

  acb_clear(I2);
  acb_clear(I6);
}
