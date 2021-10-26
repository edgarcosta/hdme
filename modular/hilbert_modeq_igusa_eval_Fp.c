
#include "modular.h"

int hilbert_modeq_igusa_eval_Fp(fmpz_mod_poly_struct* pol_vec,
				const fmpz* rs, slong ell, slong delta,
				const fmpz_mod_ctx_t ctx)
{
  fmpq* rs_lift;
  int success;
  fmpz_poly_struct num_vec[3];
  fmpz_t den;
  fmpz_t one;
  int v = MODEQ_VERBOSE;
  slong k;

  rs_lift = _fmpq_vec_init(2);
  for (k = 0; k < 3; k++) fmpz_poly_init(&num_vec[k]);
  fmpz_init(den);
  fmpz_init(one);
  
  fmpz_one(one);
  for (k = 0; k < 2; k++) fmpq_set_fmpz_frac(&rs_lift[k], &rs[k], one);

  success = hilbert_modeq_igusa_eval_Q(num_vec, den, rs_lift, ell, delta);
  if (success)
    {
      success = fmpz_mod_divides(den, one, den, ctx);
      if (v && !success)  flint_printf("(hilbert_modeq_igusa_eval_Fp) Denominator reduces to zero in the finite field\n");
    }
  if (success)
    {
      for (k = 0; k < 2; k++)
	{
	  fmpz_mod_poly_set_fmpz_poly(&pol_vec[k], &num_vec[k], ctx);
	  fmpz_mod_poly_scalar_mul_fmpz(&pol_vec[k], &pol_vec[k], den, ctx);
	}
    }

  _fmpq_vec_clear(rs_lift, 2);
  for (k = 0; k < 3; k++) fmpz_poly_clear(&num_vec[k]);
  fmpz_clear(den);
  fmpz_clear(one);

  return success;
}
