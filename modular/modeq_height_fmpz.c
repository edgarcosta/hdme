
#include "modular.h"

slong modeq_height_fmpz(const fmpz* j, slong len)
{
  slong k;
  slong h = 1;
  
  for (k = 0; k < len; k++)
    {
      if (h < fmpz_bits(&j[k])) h = fmpz_bits(&j[k]);
    }
  return h;
}
