
#include "modular.h"

int siegel_modeq_isog_igusa_Fp(fmpz* j, const fmpz_mod_poly_t pol1, const fmpz_mod_poly_t pol2,
			       const fmpz_mod_poly_t pol3, const fmpz_t root,
			       const fmpz_mod_ctx_t ctx)
{
  fmpz_mod_poly_t pol1_der;
  fmpz_t num, den;
  int res;

  fmpz_mod_poly_init(pol1_der, ctx);
  fmpz_init(num);
  fmpz_init(den);
  
  fmpz_mod_poly_derivative(pol1_der, pol1, ctx);
  fmpz_mod_poly_evaluate_fmpz(den, pol1_der, root, ctx);

  if (fmpz_is_zero(den)) res = 0;
  else
    {
      res = 1;
      fmpz_mod_inv(den, den, ctx);
      
      fmpz_set(&j[0], root);
      fmpz_mod_poly_evaluate_fmpz(num, pol2, root, ctx);
      fmpz_mod_mul(&j[1], num, den, ctx);
      fmpz_mod_poly_evaluate_fmpz(num, pol3, root, ctx);
      fmpz_mod_mul(&j[2], num, den, ctx);
    }

  fmpz_mod_poly_clear(pol1_der, ctx);
  fmpz_clear(num);
  fmpz_clear(den);
  return res;
}
