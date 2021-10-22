
#include "modular.h"

void hilbert_modeq_gundlach_simplify(fmpz_poly_t num1, fmpz_poly_t num2,
				     fmpz_t den, slong ell, slong delta)
{
  fmpz_t gcd;
  fmpz_t coeff;
  slong gcd_bits;
  slong d = 2*hilbert_nb_cosets(ell);
  slong k;

  fmpz_init(gcd);
  fmpz_init(coeff);
  
  fmpz_set(gcd, den);
  for (k = 0; k < d + 1; k++)
    {
      fmpz_poly_get_coeff_fmpz(coeff, num1, k);
      fmpz_gcd(gcd, gcd, coeff);
      fmpz_poly_get_coeff_fmpz(coeff, num2, k);
      fmpz_gcd(gcd, gcd, coeff);
    }
  
  gcd_bits = fmpz_bits(gcd);
  flint_printf("(hilbert_modeq_gundlach_simplify) Simplify result by %wd-bit integer\n",
	       gcd_bits);  
  fmpz_divexact(den, den, gcd);
  fmpz_poly_scalar_divexact_fmpz(num1, num1, gcd);
  fmpz_poly_scalar_divexact_fmpz(num2, num2, gcd);

  fmpz_clear(gcd);
  fmpz_clear(coeff);
}
