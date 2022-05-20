
#include "modular.h"

int siegel_modeq_eval_Fp(fmpz_mod_poly_struct* pol_vec,
			 const fmpz* j, slong ell, const fmpz_mod_ctx_t ctx)
{
  fmpq* j_lift;
  fmpz_poly_struct num_vec[3];
  fmpz_t den;
  slong k;
  int success;

  j_lift = _fmpq_vec_init(3);
  for (k = 0; k < 3; k++) fmpz_poly_init(&num_vec[k]);
  fmpz_init(den);

  modeq_input_lift(j_lift, j, 3);
  success = siegel_modeq_eval_Q(num_vec, den, j_lift, ell);
  if (success) success = modeq_reduce(pol_vec, num_vec, den, 3, ctx);

  _fmpq_vec_clear(j_lift, 3);
  for (k = 0; k < 3; k++) fmpz_poly_clear(&num_vec[k]);
  fmpz_clear(den);
  return success;
}
