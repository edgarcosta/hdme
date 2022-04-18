
#include "igusa.h"

void cov_divexact_fmpz_si(fmpz* I, fmpz* S, slong scal)
{
  fmpz_t f;
  fmpz_init(f);

  fmpz_set_si(f, scal);
  cov_divexact_fmpz(I, f, S);
  
  fmpz_clear(f);
}
