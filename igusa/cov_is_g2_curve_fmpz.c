
#include "igusa.h"

int cov_is_g2_curve_fmpz(fmpz* I)
{
  return !fmpz_is_zero(cov_I10(I));
}
