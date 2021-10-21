
#include "hilbert.h"


int hilbert_modeq_sym_igusa_Q(fmpz_poly_t num1, fmpz_poly_t num2, fmpz_poly_t num3,
			      fmpz_t den, const acb_poly_t pol1_acb,
			      const acb_poly_t pol2_acb, const acb_poly_t pol3_acb,
			      slong ell, slong prec)
{
  fmpq_t c;
  fmpz_t c_num;
  fmpz_t current_den;
  fmpq_poly_t pol1, pol2, pol3;
  acb_t x;
  slong n = hilbert_nb_cosets(ell, delta);
  int success = 1;
  slong k;

  fmpq_init(c);
  fmpz_init(c_num);
  fmpz_init(current_den);
  fmpq_poly_init(pol1);
  fmpq_poly_init(pol2);
  fmpq_poly_init(pol3);
  acb_init(x);

  fmpz_one(current_den);

  for (k = 0; (k <= 2*n) && success; k++)
    {
      acb_poly_get_coeff_acb(x, pol1_acb, k);
      success = hilbert_modeq_coeff_Q(c, current_den, x, current_den, prec);
      if (success)
	{
	  fmpq_poly_set_coeff_fmpq(pol1, k, c);
	}
    }
  
  for (k = 0; (k <= 2*n-1) && success; k++)
    {
      acb_poly_get_coeff_acb(x, pol2_acb, k);
      success = hilbert_modeq_coeff_Q(c, current_den, x, current_den, prec);
      if (success)
	{
	  fmpq_poly_set_coeff_fmpq(pol2, k, c);
	}
    }
  
  for (k = 0; (k <= 2*n-1) && success; k++)
    {
      acb_poly_get_coeff_acb(x, pol3_acb, k);
      success = hilbert_modeq_coeff_Q(c, current_den, x, current_den, prec);
      if (success)
	{
	  fmpq_poly_set_coeff_fmpq(pol3, k, c);
	}
    }

  if (success)
    {
      fmpz_set(den, current_den);
      fmpq_poly_scalar_mul_fmpz(pol1, pol1, den);
      fmpq_poly_scalar_mul_fmpz(pol2, pol2, den);
      fmpq_poly_scalar_mul_fmpz(pol3, pol3, den);
      for (k = 0; k <= 2*n-1; k++)
	{
	  fmpq_poly_get_coeff_fmpz(c_num, pol1, k, c);
	  fmpz_poly_set_coeff_fmpz(num1, k, c_num);
	  fmpq_poly_get_coeff_fmpz(c_num, pol2, k, c);
	  fmpz_poly_set_coeff_fmpz(num2, k, c_num);
	  fmpq_poly_get_coeff_fmpz(c_num, pol2, k, c);
	  fmpz_poly_set_coeff_fmpz(num2, k, c_num);
	}
      fmpq_poly_get_coeff_fmpz(c_num, pol1, 2*n, c);
      fmpz_poly_set_coeff_fmpz(num1, 2*n, c_num);
    }

  fmpq_clear(c);
  fmpz_clear(c_num);
  fmpz_clear(current_den);
  fmpq_poly_clear(pol1);
  fmpq_poly_clear(pol2);
  fmpq_poly_clear(pol3);
  acb_clear(x);

  return success;
}
