
#include "igusa.h"

void igusa_switch_I6_fmpz(fmpz* I, fmpz* S)
{
  int r;

  r = igusa_I6_fmpz(&I[2], S);
  if (!r)
    {
      cov_rescale_fmpz_si(I, S, 3);
      igusa_I6_fmpz(&I[2], I);
    }
  cov_normalize_fmpz(I, I);
}
