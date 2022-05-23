
#include "modular.h"

void modeq_isog_monomials_Fp(fmpz* M, const modeq_t E,
			     const fmpz_t root, slong mult, const fmpz_mod_ctx_t ctx)
{
  fmpz_mod_poly_t num;
  
  fmpz_mod_poly_init(num, ctx);
  
  for (j = 0; k < modeq_nb(E); k++)
    {
      fmpz_mod_poly_set_fmpz_poly(num, modeq_interpolate(E, j), ctx);
      pol_remove_root_Fp(num, num, root, mult-1);
      fmpz_mod_poly_evaluate_fmpz(&M[j], num, root);
    }

  fmpz_mod_poly_clear(num, ctx);
}
