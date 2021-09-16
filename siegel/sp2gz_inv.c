
#include "siegel.h"

void
sp2gz_inv(sp2gz_t r, const sp2gz_t m)
{
  if (r != m)
    sp2gz_set(r, m);
  
  fmpz_mat_transpose(&r->a, &r->a);
  fmpz_mat_transpose(&r->b, &r->b);
  fmpz_mat_transpose(&r->c, &r->c);
  fmpz_mat_transpose(&r->d, &r->d);
  fmpz_mat_swap(&r->a, &r->d);
  fmpz_mat_neg(&r->b, &r->b);
  fmpz_mat_neg(&r->c, &r->c);
}
