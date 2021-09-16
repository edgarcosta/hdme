
#include "siegel.h"

void
siegel_star(acb_mat_t w, const sp2gz_t m, const acb_mat_t z, slong prec)
{
  acb_mat_t cd, x;
  acb_mat_init(x, m->g, m->g);
  acb_mat_init(cd, m->g, m->g);

  acb_mat_set_fmpz_mat(cd, &m->c);
  acb_mat_mul(x, cd, z, prec);
  acb_mat_set_fmpz_mat(cd, &m->d);
  acb_mat_add(x, x, cd, prec);

  acb_mat_set(w, x);

  acb_mat_clear(cd);
  acb_mat_clear(x);
}
