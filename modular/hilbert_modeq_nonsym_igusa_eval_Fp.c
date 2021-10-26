
#include "modular.h"


int hilbert_modeq_nonsym_igusa_eval_Fp(fmpz_mod_poly_struct* pol_vec,
				       const fmpz* rs, slong ell, const fmpz_poly_t beta,
				       slong delta, const fmpz_mod_ctx_t ctx)
{
  fmpq* rs_lift;
  int success;
  fmpz_poly_struct num_vec[3];
  fmpz_t den;
  slong k;

  rs_lift = _fmpq_vec_init(2);
  for (k = 0; k < 3; k++) fmpz_poly_init(&num_vec[k]);
  fmpz_init(den);

  modeq_input_lift(rs_lift, rs, 2);
  success = hilbert_modeq_nonsym_igusa_eval_Q(num_vec, den, rs_lift, ell, beta, delta);
  if (success) success = modeq_reduce(pol_vec, num_vec, den, 3, ctx);
  
  _fmpq_vec_clear(rs_lift, 2);
  for (k = 0; k < 3; k++) fmpz_poly_clear(&num_vec[k]);
  fmpz_clear(den);

  return success;
}
