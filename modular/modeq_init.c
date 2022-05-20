
#include "modular.h"

void modeq_init(modeq_t E, slong nb)
{
  slong k;
  modeq_nb(E) = nb;
  fmpz_poly_init(modeq_equation(E));
  fmpz_init(modeq_den(E));
  E->interp = flint_malloc(nb * sizeof(fmpz_poly_struct));
  for (k = 0; k < nb; k++) fmpz_poly_init(modeq_interpolate(E, k));
}
