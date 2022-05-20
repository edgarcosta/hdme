
#include "igusa.h"

int igusa_is_g2_curve(acb_srcptr I)
{
  return !acb_contains_zero(igusa_I10(I));
}
