
#include "modular.h"

void modeq_simplify(fmpz_poly_struct* num_vec, fmpz_t den, slong degree, slong nb)
{
  fmpz_t gcd;
  fmpz_t coeff;
  slong gcd_bits;
  slong i, k;
  int v = MODEQ_VERBOSE;
  
  fmpz_init(gcd);
  fmpz_init(coeff);
  
  fmpz_set(gcd, den);
  for (i = 0; i < nb; i++)
    {
      for (k = 0; k <= degree ; k++)
	{
	  fmpz_poly_get_coeff_fmpz(coeff, &num_vec[i], k);
	  fmpz_gcd(gcd, gcd, coeff);
	}
    }
  
  fmpz_divexact(den, den, gcd);
  for (i = 0; i < nb; i++)
    {
      fmpz_poly_scalar_divexact_fmpz(&num_vec[i], &num_vec[i], gcd);
    }
  
  gcd_bits = fmpz_bits(gcd);
  if (v) flint_printf("(modeq_simplify) Simplify result by %wd-bit integer\n", gcd_bits);
  
  fmpz_clear(gcd);
  fmpz_clear(coeff);
}

void siegel_modeq_simplify(fmpz_poly_t num1, fmpz_poly_t num2,
			   fmpz_poly_t num3, fmpz_t den, slong ell)
{
  fmpz_t gcd;
  fmpz_t coeff;
  slong gcd_bits;
  slong k;

  fmpz_init(gcd);
  fmpz_init(coeff);
  
  fmpz_set(gcd, den);
  for (k = 0; k < siegel_nb_cosets(ell) + 1; k++)
    {
      fmpz_poly_get_coeff_fmpz(coeff, num1, k);
      fmpz_gcd(gcd, gcd, coeff);
      fmpz_poly_get_coeff_fmpz(coeff, num2, k);
      fmpz_gcd(gcd, gcd, coeff);
      fmpz_poly_get_coeff_fmpz(coeff, num3, k);
      fmpz_gcd(gcd, gcd, coeff);
    }
  
  gcd_bits = fmpz_bits(gcd);
  flint_printf("(siegel_modeq_simplify) Simplify result by %wd-bit integer\n", gcd_bits);
  /*
  fmpz_print(gcd);
  flint_printf("\n"); */
  
  fmpz_divexact(den, den, gcd);
  fmpz_poly_scalar_divexact_fmpz(num1, num1, gcd);
  fmpz_poly_scalar_divexact_fmpz(num2, num2, gcd);
  fmpz_poly_scalar_divexact_fmpz(num3, num3, gcd);

  fmpz_clear(gcd);
  fmpz_clear(coeff);
}
