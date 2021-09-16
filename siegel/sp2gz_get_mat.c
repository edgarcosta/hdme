
#include "siegel.h"

void
sp2gz_get_mat(fmpz_mat_t m, const sp2gz_t n)
{
  fmpz_mat_t ac, bd;
  slong g = n->g;

  fmpz_mat_init(ac, 2*g, g);
  fmpz_mat_init(bd, 2*g, g);

  fmpz_mat_concat_vertical(ac, &n->a, &n->c);
  fmpz_mat_concat_vertical(bd, &n->b, &n->d);
  fmpz_mat_concat_horizontal(m, ac, bd);

  fmpz_mat_clear(ac);
  fmpz_mat_clear(bd);
}
