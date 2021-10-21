
#include "modular.h"

slong hilbert_modeq_height(fmpq* params, slong len)
{
  fmpz_t den;
  fmpz* num;
  slong k;
  slong h;

  fmpz_init(den);
  num = _fmpz_vec_init(len);
  
  siegel_modeq_fmpz_input(den, num, params, len);
  h = fmpz_bits(den);
  for (k = 0; k < len; k++)
    {
      if (h < fmpz_bits(&num[k])) h = fmpz_bits(&num[k]);
    }
  
  fmpz_clear(den);
  _fmpz_vec_clear(num, len);
  return h;
}
