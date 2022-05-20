
#include "igusa.h"

int igusa_I2_fmpz(fmpz_t I2, fmpz* I)
{
  if (!fmpz_divisible(igusa_I12(I), igusa_I10(I)))
    {
      return 0;
    }
  else
    {
      fmpz_divexact(I2, igusa_I12(I), igusa_I10(I));
      return 1;
    }      
}
