
#include "hilbert.h"

void hilbert_sigma1(acb_t z, const fmpz_poly_t x, slong delta, slong prec)
{
  acb_t res;
  acb_init(res);

  acb_set_fmpz(z, fmpz_poly_get_coeff_ptr(x, 0));

  acb_zero(res);
  arb_sqrt_ui(acb_realref(res), delta, prec);
  if (delta % 2 == 1)
    {
      acb_add_si(res, res, 1, prec);
    }
  acb_div_si(res, res, 2, prec);
  acb_mul_fmpz(res, res, fmpz_poly_get_coeff_ptr(x, 1), prec);
	       
  acb_add(z, z, res, prec);
  
  acb_clear(res);
}
