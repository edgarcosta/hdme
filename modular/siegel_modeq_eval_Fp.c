
#include "modular.h"

int siegel_modeq_eval_Fp(fmpz_mod_poly_struct* pol_vec,
			 const fmpz* j, slong ell, const fmpz_mod_ctx_t ctx)
{
  fmpq* j_lift;
  fmpz_poly_struct num_vec[3];
  fmpz_t den;
  fmpz_t one;
  slong k;
  int success;
  int v = MODEQ_VERBOSE;

  j_lift = _fmpq_vec_init(3);
  for (k = 0; k < 3; k++) fmpz_poly_init(&num_vec[k]);
  fmpz_init(den);
  fmpz_init(one);

  fmpz_one(one);
  for (k = 0; k < 3; k++) fmpq_set_fmpz_frac(&j_lift[k], &j[k], one);

  success = siegel_modeq_eval_Q(num_vec, den, j_lift, ell);
  if (success)
    {
      success = fmpz_mod_divides(den, one, den, ctx);
      if (v && !success) flint_printf("(siegel_modeq_eval_Fp) Denominator vanishes\n");
    }
  if (success)
    {
      for (k = 0; k < 3; k++)
	{
	  fmpz_mod_poly_set_fmpz_poly(&pol_vec[k], &num_vec[k], ctx);
	  fmpz_mod_poly_scalar_mul_fmpz(&pol_vec[k], &pol_vec[k], den, ctx);
	}
    }

  _fmpq_vec_clear(j_lift, 3);
  for (k = 0; k < 3; k++) fmpz_poly_clear(&num_vec[k]);
  fmpz_clear(den);
  fmpz_clear(one);
  return success;
}
