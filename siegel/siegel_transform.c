
#include "siegel.h"

int
siegel_transform(acb_mat_t w, const sp2gz_t m, const acb_mat_t z, slong prec)
{
  acb_mat_t abcd, num, den, invden;
  int valid_output;
  
  acb_mat_init(abcd, m->g, m->g);
  acb_mat_init(num, m->g, m->g);
  acb_mat_init(den, m->g, m->g);
  acb_mat_init(invden, m->g, m->g);

  acb_mat_set_fmpz_mat(abcd, &m->a);
  acb_mat_mul(num, abcd, z, prec);
  acb_mat_set_fmpz_mat(abcd, &m->b);
  acb_mat_add(num, num, abcd, prec);

  acb_mat_set_fmpz_mat(abcd, &m->c);
  acb_mat_mul(den, abcd, z, prec);
  acb_mat_set_fmpz_mat(abcd, &m->d);
  acb_mat_add(den, den, abcd, prec);

  valid_output = acb_mat_inv(invden, den, prec);

  if (!valid_output)
    {
      flint_fprintf(stderr, "Impossible inverse:\n");
      acb_mat_fprintd(stderr, den, 30);
      flint_fprintf(stderr, "\n\n");
    }
    
  acb_mat_mul(w, num, invden, prec);
  
  acb_mat_clear(abcd);
  acb_mat_clear(num);
  acb_mat_clear(den);
  acb_mat_clear(invden);
  return valid_output;
}
