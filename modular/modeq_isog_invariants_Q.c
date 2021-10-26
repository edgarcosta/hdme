
#include "modular.h"

int modeq_isog_invariants_Q(fmpq* j, const fmpz_poly_struct* num_vec,
			    const fmpq_t root, slong nb)
{
  fmpz_poly_t num1_der;
  fmpq_t num, den;
  int res;
  slong k;
  
  fmpz_poly_init(num1_der);
  fmpq_init(num);
  fmpq_init(den);
  
  fmpz_poly_derivative(num1_der, &num_vec[0]);
  fmpz_poly_evaluate_fmpq(den, num1_der, root);
  
  if (fmpq_is_zero(den)) res = 0;
  else
    {
      res = 1;
      fmpq_set(&j[0], root);
      for (k = 1; k < nb; k++)
	{
	  fmpz_poly_evaluate_fmpq(num, &num_vec[k], root);
	  fmpq_div(&j[k], num, den);
	}
    }

  fmpz_poly_clear(num1_der);
  fmpq_clear(num);
  fmpq_clear(den);  
  return res;  
}
