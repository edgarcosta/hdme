
#include "igusa.h"

int igusa_I2_fmpz(fmpz_t I2, fmpz* I)
{
  if (!fmpz_divisible(cov_I12(I), cov_I10(I)))
    {
      return 0;
    }
  else
    {
      fmpz_divexact(I2, cov_I12(I), cov_I10(I));
      return 1;
    }      
}
