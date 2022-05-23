
#include "modular.h"

int hilbert_isog_monomials_Fp(slong* nb_roots, fmpz* all_M, slong nb_M,
			      slong* exp_array, fmpz* I, slong ell,
			      slong delta, const fmpz_mod_ctx_t ctx)
{
  modeq_t E;
  modeq_ctx_t ctx;
  int res;

  modeq_init(E);
  modeq_ctx_init(ctx);
  
  res = hilbert_modeq_eval(E, ctx, I, ell, delta);
  if (res) res = modeq_all_isog_Fp(nb_roots, all_M, nb_M, exp_array,
				   E, ctx, fpctx);
  
  modeq_clear(E);
  modeq_ctx_clear(ctx);
  return res;
}
