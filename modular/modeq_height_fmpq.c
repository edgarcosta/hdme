
#include "modular.h"

slong modeq_height_fmpq(fmpq* j, slong len)
{
  fmpz_t den;
  fmpz* num;
  slong k;
  slong h;

  fmpz_init(den);
  num = _fmpz_vec_init(len);

  modeq_input_get_fmpz(den, num, j, len);
  h = fmpz_bits(den);
  for (k = 0; k < len; k++)
    {
      if (h < fmpz_bits(&num[k])) h = fmpz_bits(&num[k]);
    }

  fmpz_clear(den);
  _fmpz_vec_clear(num, len);
  return h;
}
