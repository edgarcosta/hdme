
#include "modular.h"

void modeq_input_lift(fmpq* j, const fmpz* input, slong nb)
{
  fmpz_t one;
  slong k;
  
  fmpz_init(one);
  
  fmpz_one(one);
  for (k = 0; k < nb; k++)
    {
      fmpq_set_fmpz_frac(&j[k], &input[k], one);
    }

  fmpz_clear(one);
}
