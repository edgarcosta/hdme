
#include "modular.h"

int siegel_modeq_isog_igusa_Q(fmpq* j, const fmpz_poly_t num1, const fmpz_poly_t num2,
			      const fmpz_poly_t num3, const fmpq_t root)
{
  fmpz_poly_t num1_der;
  fmpq_t num, den;
  int res;

  fmpz_poly_init(num1_der);
  fmpq_init(num);
  fmpq_init(den);
  
  fmpz_poly_derivative(num1_der, num1);
  fmpz_poly_evaluate_fmpq(den, num1_der, root);

  if (fmpq_is_zero(den)) res = 0;
  else
    {
      res = 1;
      fmpq_set(&j[0], root);
      fmpz_poly_evaluate_fmpq(num, num2, root);
      fmpq_div(&j[1], num, den);
      fmpz_poly_evaluate_fmpq(num, num3, root);
      fmpq_div(&j[2], num, den);
    }

  fmpz_poly_clear(num1_der);
  fmpq_clear(num);
  fmpq_clear(den);
  return res;
}
