
#include "modular.h"

int hilbert_modeq_gundlach_eval_Fp(fmpz_mod_poly_t pol1, fmpz_mod_poly_t pol2,
				   const fmpz* g, slong ell,
				   slong delta, const fmpz_mod_ctx_t ctx)
{  
  fmpq* g_lift;
  int success;
  fmpz_poly_t num1, num2;
  fmpz_t den;
  fmpz_t one;
  int v = MODEQ_VERBOSE;
  
  g_lift = _fmpq_vec_init(2);
  fmpz_poly_init(num1);
  fmpz_poly_init(num2);
  fmpz_init(den);
  fmpz_init(one);
  
  fmpz_one(one);
  fmpq_set_fmpz_frac(&g_lift[0], &g[0], one);
  fmpq_set_fmpz_frac(&g_lift[1], &g[1], one);
  
  success = hilbert_modeq_gundlach_eval_Q(num1, num2, den, g_lift, ell, delta);
  if (success)
    {
      success = fmpz_mod_divides(den, one, den, ctx);
      if (v && !success)  flint_printf("(hilbert_modeq_gundlach_eval_Fp) Denominator reduces to zero in the finite field\n");
    }
  if (success)
    {
      fmpz_mod_poly_set_fmpz_poly(pol1, num1, ctx);
      fmpz_mod_poly_scalar_mul_fmpz(pol1, pol1, den, ctx);
      fmpz_mod_poly_set_fmpz_poly(pol2, num2, ctx);
      fmpz_mod_poly_scalar_mul_fmpz(pol2, pol2, den, ctx);
    }
  
  _fmpq_vec_clear(g_lift, 2);
  fmpz_poly_clear(num1);
  fmpz_poly_clear(num2);
  fmpz_clear(den);
  fmpz_clear(one);

  return success;
}
