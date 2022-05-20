
#include "modular.h"

void modeq_acb_init(modeq_acb_t E, slong nb)
{
  slong k;
  modeq_nb(E) = nb;
  acb_poly_init(modeq_equation(E));
  acb_init(modeq_den(E));
  E->interp = flint_malloc(nb * sizeof(acb_poly_struct));
  for (k = 0; k < nb; k++) acb_poly_init(modeq_interpolate(E, k));
}
