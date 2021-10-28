
#include "siegel.h"

int siegel_transform(acb_mat_t w, const fmpz_mat_t m, const acb_mat_t z, slong prec)
{
  slong g = fmpz_mat_half_dim(m);
  fmpz_mat_t a;
  acb_mat_t x, num, den, invden;
  int valid_output;

  fmpz_mat_init(a, g, g);  
  acb_mat_init(x, g, g);
  acb_mat_init(num, g, g);
  acb_mat_init(den, g, g);
  acb_mat_init(invden, g, g);

  fmpz_mat_get_a(a, m);
  acb_mat_set_fmpz_mat(x, a);
  acb_mat_mul(num, x, z, prec);
  fmpz_mat_get_b(a, m);
  acb_mat_set_fmpz_mat(x, a);
  acb_mat_add(num, num, x, prec);

  siegel_star(den, m, z, prec);
  valid_output = acb_mat_inv(invden, den, prec);
  if (!valid_output)
    {
      flint_fprintf(stderr, "(siegel_transform) Warning: impossible inverse\n");
      acb_mat_fprintd(stderr, den, 30);
      flint_fprintf(stderr, "\n");
    }
    
  acb_mat_mul(w, num, invden, prec);

  fmpz_mat_clear(a);
  acb_mat_clear(x);
  acb_mat_clear(num);
  acb_mat_clear(den);
  acb_mat_clear(invden);
  return valid_output;
}
