
#include "modular.h"

int hilbert_modeq_gundlach_eval_Fp(fmpz_mod_poly_struct* pol_vec,
				   const fmpz* g, slong ell,
				   slong delta, const fmpz_mod_ctx_t ctx)
{  
  fmpq* g_lift;
  int success;
  fmpz_poly_struct num_vec[2];
  fmpz_t den;
  slong k;
  
  g_lift = _fmpq_vec_init(2);
  for (k = 0; k < 2; k++) fmpz_poly_init(&num_vec[k]);
  fmpz_init(den);

  modeq_input_lift(g_lift, g, 2);  
  success = hilbert_modeq_gundlach_eval_Q(num_vec, den, g_lift, ell, delta);
  if (success) success = modeq_reduce(pol_vec, num_vec, den, 2, ctx);
  
  _fmpq_vec_clear(g_lift, 2);
  for (k = 0; k < 2; k++) fmpz_poly_clear(&num_vec[k]);
  fmpz_clear(den);

  return success;
}
