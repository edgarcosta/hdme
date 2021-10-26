
#include "modular.h"

int siegel_modeq_eval_Fp(fmpz_mod_poly_t pol1, fmpz_mod_poly_t pol2, fmpz_mod_poly_t pol3,
			 const fmpz* j, slong ell, const fmpz_mod_ctx_t ctx)
{
  fmpq* j_lift;
  fmpz_poly_t num1, num2, num3;
  fmpz_t den;
  fmpz_t one;
  slong k;
  int success;
  int v = MODEQ_VERBOSE;

  j_lift = _fmpq_vec_init(3);
  fmpz_poly_init(num1);
  fmpz_poly_init(num2);
  fmpz_poly_init(num3);
  fmpz_init(den);
  fmpz_init(one);

  fmpz_one(one);
  for (k = 0; k < 2; k++) fmpq_set_fmpz_frac(&j_lift[k], &j[k], one);

  success = siegel_modeq_eval_Q(num1, num2, num3, den, j_lift, ell);
  if (success)
    {
      success = fmpz_mod_divides(den, one, den, ctx);
      if (v && !success) flint_printf("(siegel_modeq_eval_Fp) Denominator vanishes\n");
    }
  if (success)
    {
      fmpz_mod_poly_set_fmpz_poly(pol1, num1, ctx);
      fmpz_mod_poly_set_fmpz_poly(pol2, num2, ctx);
      fmpz_mod_poly_set_fmpz_poly(pol3, num3, ctx);
      fmpz_mod_poly_scalar_mul_fmpz(pol1, pol1, den, ctx);
      fmpz_mod_poly_scalar_mul_fmpz(pol2, pol2, den, ctx);
      fmpz_mod_poly_scalar_mul_fmpz(pol3, pol3, den, ctx);
    }

  _fmpq_vec_clear(j_lift, 3);
  fmpz_poly_clear(num1);
  fmpz_poly_clear(num2);
  fmpz_poly_clear(num3);
  fmpz_clear(den);
  fmpz_clear(one);
  return success;
}
