
#include "modular.h"

slong siegel_modeq_height_fmpz(const fmpz* j)
{
  slong k;
  slong h = 1;
  fmpz* j_mod;

  j_mod = _fmpz_vec_init(3);

  _fmpz_vec_set(j_mod, j, 3);
  fmpz_mul(&j_mod[0], &j_mod[0], &j_mod[0]);
  
  for (k = 0; k < 3; k++)
    {
      if (h < fmpz_bits(&j[k])) h = fmpz_bits(&j[k]);
    }
  return h;

  _fmpz_vec_clear(j_mod, 3);
}
