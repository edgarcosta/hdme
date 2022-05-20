
#include "modular.h"

void modeq_clear(modeq_t E)
{
  slong k;
  fmpz_poly_clear(modeq_equation(E));
  fmpz_clear(modeq_den(E));
  for (k = 0; k < nb; k++) fmpz_poly_clear(modeq_interpolate(E, k));
  flint_free(E->interp);
}
