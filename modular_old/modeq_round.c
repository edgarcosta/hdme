
#include "modular.h"

int modeq_round(fmpz_poly_struct* num_vec, fmpz_t den, const acb_poly_struct* num_acb_vec,
		const acb_t den_acb, slong degree, slong nb)
{  
  acb_t coeff;
  arf_t radius, max_radius;
  fmpz_t rd;
  slong radius_prec = 100;
  slong radius_bits;
  int res = 1;
  slong k;
  slong v = MODEQ_VERBOSE;
  
  acb_init(coeff);
  fmpz_init(rd);
  arf_init(radius);
  arf_init(max_radius);
  
  arf_zero(max_radius);
  acb_get_rad_ubound_arf(radius, den_acb, radius_prec);
  arf_max(max_radius, max_radius, radius);  
  res = modeq_round_coeff(den, den_acb);

  for (k = 0; k < nb; k++)
    {
      if (res)
	{
	  res = modeq_round_poly(&num_vec[k], radius,
				 &num_acb_vec[k], degree);
	  arf_max(max_radius, max_radius, radius);	  
	}
    }

  radius_bits = arf_abs_bound_lt_2exp_si(max_radius);
  if (v && res) flint_printf("(modeq_round) Excess precision: %wd\n", -radius_bits);
  
  acb_clear(coeff);
  fmpz_clear(rd);
  arf_clear(radius);
  arf_clear(max_radius);
  return res;  
}
