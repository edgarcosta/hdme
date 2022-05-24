
#include "igusa.h"

void cov_normalize_fmpz_wt1(fmpz* I, fmpz* S, slong nb)
{
  fmpz_t g;

  fmpz_init(g);
  _fmpz_vec_content(g, S, nb);
  _fmpz_vec_scalar_divexact_fmpz(I, S, nb, g);
  fmpz_clear(g);
}
