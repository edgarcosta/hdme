
#include "modular.h"

int hilbert_modeq_sym_igusa_eval_Fp(fmpz_mod_poly_t pol1, fmpz_mod_poly_t pol2,
				    fmpz_mod_poly_t pol3,
				    const fmpz* rs, slong ell, slong delta,
				    const fmpz_mod_ctx_t ctx)
{
  fmpq* rs_lift;
  int success;
  fmpz_poly_t num1, num2, num3;
  fmpz_t den;
  fmpz_t one;
  int v = HILBERT_VERBOSE;

  rs_lift = _fmpq_vec_init(2);
  fmpz_poly_init(num1);
  fmpz_poly_init(num2);
  fmpz_poly_init(num3);
  fmpz_init(den);
  fmpz_init(one);
  
  fmpz_one(one);
  fmpq_set_fmpz_frac(&rs_lift[0], &rs[0], one);
  fmpq_set_fmpz_frac(&rs_lift[1], &rs[1], one);

  success = hilbert_modeq_sym_igusa_eval_Q(num1, num2, num3, den, rs_lift, ell, delta);
  if (success)
    {
      success = fmpz_mod_divides(den, one, den, ctx);
      if (v && !success)  flint_printf("(hilbert_modeq_sym_igusa_eval_Fp) Denominator reduces to zero in the finite field\n");
    }
  if (success)
    {
      fmpz_mod_poly_set_fmpz_poly(pol1, num1, ctx);
      fmpz_mod_poly_scalar_mul_fmpz(pol1, pol1, den, ctx);
      fmpz_mod_poly_set_fmpz_poly(pol2, num2, ctx);
      fmpz_mod_poly_scalar_mul_fmpz(pol2, pol2, den, ctx);
      fmpz_mod_poly_set_fmpz_poly(pol3, num3, ctx);
      fmpz_mod_poly_scalar_mul_fmpz(pol3, pol3, den, ctx);
    }

  _fmpq_vec_clear(rs_lift, 2);
  fmpz_poly_clear(num1);
  fmpz_poly_clear(num2);
  fmpz_poly_clear(num3);
  fmpz_clear(den);
  fmpz_clear(one);

  return success;
}
