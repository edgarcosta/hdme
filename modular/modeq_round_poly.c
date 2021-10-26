
#include "modular.h"

int modeq_round_poly(fmpz_poly_t pol, arf_t max_radius,
		     const acb_poly_t pol_acb, slong degree)
{
  
  acb_t coeff;
  arf_t radius;
  fmpz_t rd;
  slong k;
  slong radius_prec = 100;
  int res = 1;

  acb_init(coeff);
  fmpz_init(rd);
  arf_init(radius);
  
  arf_zero(max_radius);
  fmpz_poly_zero(pol);
  
  for (k = 0; k < degree + 1; k++)
    {
      fmpz_zero(rd);
      if (res)
	{
	  acb_poly_get_coeff_acb(coeff, pol_acb, k);
	  res = modeq_round_coeff(rd, coeff);
	  acb_get_rad_ubound_arf(radius, coeff, radius_prec);
	  arf_max(max_radius, max_radius, radius);
	  fmpz_poly_set_coeff_fmpz(pol, k, rd);
	}
    }

  acb_clear(coeff);
  fmpz_clear(rd);
  arf_clear(radius);
  return res;
}
