
#include "modular.h"

int modeq_isog_invariants_C(acb_ptr j, const fmpz_poly_struct* num_vec,
			    const acb_t root, slong nb, slong prec)
{
  acb_poly_struct nums[3];
  acb_poly_t num1_der;
  acb_t num, den;
  int res;
  slong k;

  for (k = 0; k < 3; k++) acb_poly_init(&nums[k]);
  acb_poly_init(num1_der);
  acb_init(num);
  acb_init(den);

  for (k = 0; k < 3; k++) acb_poly_set_fmpz_poly(&nums[k], &num_vec[k], prec);
  acb_poly_derivative(num1_der, &nums[0], prec);
  acb_poly_evaluate(den, num1_der, root, prec);
  
  if (acb_contains_zero(den)) res = 0;
  else
    {
      res = 1;
      acb_set(&j[0], root);
      for (k = 1; k < nb; k++)
	{
	  acb_poly_evaluate(num, &nums[k], &j[0], prec);
	  acb_div(&j[k], num, den, prec);
	}
    }

  for (k = 0; k < 3; k++) acb_poly_clear(&nums[k]);
  acb_poly_clear(num1_der);
  acb_clear(num);
  acb_clear(den);  
  return res;  
}
