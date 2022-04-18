
#include "igusa.h"

void cov_rescale_fmpz_si(fmpz* I, fmpz* S, slong scal)
{
  fmpz_t f;

  fmpz_init(f);
  fmpz_set_si(f, scal);
  cov_rescale_fmpz(I, S, f);
  
  fmpz_clear(f);
}
