
#include "modular.h"

int modeq_isog_invariants_Fp(fmpz* j, const fmpz_mod_poly_struct* num_vec,
			     const fmpz_t root, slong nb,
			     const fmpz_mod_ctx_t ctx)
{
  fmpz_mod_poly_t pol1_der;
  fmpz_t num, den;
  int res;
  slong k;
  
  fmpz_mod_poly_init(pol1_der, ctx);
  fmpz_init(num);
  fmpz_init(den);
  
  fmpz_mod_poly_derivative(pol1_der, &num_vec[0], ctx);
  fmpz_mod_poly_evaluate_fmpz(den, pol1_der, root, ctx);

  if (fmpz_is_zero(den)) res = 0;
  else
    {
      res = 1;
      fmpz_mod_inv(den, den, ctx);      
      fmpz_set(&j[0], root);
      for (k = 1; k < nb; k++)
	{
	  fmpz_mod_poly_evaluate_fmpz(num, &num_vec[k], root, ctx);
	  fmpz_mod_mul(&j[k], num, den, ctx);
	}
    }
  
  fmpz_mod_poly_clear(pol1_der, ctx);
  fmpz_clear(num);
  fmpz_clear(den);
  return res;  
}
