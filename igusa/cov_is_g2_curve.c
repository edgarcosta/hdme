
#include "igusa.h"

int cov_is_g2_curve(acb_srcptr I)
{
  return !acb_contains_zero(cov_I10(I));
}
