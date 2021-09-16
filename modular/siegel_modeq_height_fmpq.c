
#include "modular.h"

slong siegel_modeq_height_fmpq(fmpq* j)
{
  fmpz_t den;
  fmpz* num;
  fmpq* j_mod;
  slong k;
  slong h;

  fmpz_init(den);
  num = _fmpz_vec_init(3);
  j_mod = _fmpq_vec_init(3);
  
  fmpq_mul(&j_mod[0], &j[0], &j[0]);
  fmpq_set(&j_mod[1], &j[1]);
  fmpq_set(&j_mod[2], &j[2]);

  siegel_modeq_fmpz_input(den, num, j_mod, 3);
  h = fmpz_bits(den);
  for (k = 0; k < 3; k++)
    {
      if (h < fmpz_bits(&num[k])) h = fmpz_bits(&num[k]);
    }

  fmpz_clear(den);
  _fmpz_vec_clear(num, 3);
  _fmpq_vec_clear(j_mod, 3);
  return h;
}
