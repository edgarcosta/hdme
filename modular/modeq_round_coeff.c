
#include "modular.h"

int modeq_round_coeff(fmpz_t c, const acb_t x)
{
  fmpz_t re, im;
  int res;
  
  fmpz_init(re);
  fmpz_init(im);

  res = arb_get_unique_fmpz(re, acb_realref(x))
    && arb_get_unique_fmpz(im, acb_imagref(x));
  if (res) res = fmpz_is_zero(im);
  if (res) fmpz_set(c, re);

  fmpz_clear(re);
  fmpz_clear(im);
  return res;
}
