
#include "modular.h"

void modeq_acb_clear(modeq_acb_t E)
{
  slong k;
  acb_poly_clear(modeq_equation(E));
  acb_clear(modeq_den(E));
  for (k = 0; k < modeq_nb(E); k++) acb_poly_clear(modeq_interpolate(E, k));
  flint_free(E->interp);
}
