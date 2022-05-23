
#include "modular.h"

int siegel_isog_monomials_Q(slong* nb_roots, fmpz* all_M, slong* nb_M,
			    slong* exp_array, fmpz* I, slong ell)
{
  modeq_t E;
  modeq_ctx_t ctx;
  int res;

  modeq_init(E);
  modeq_ctx_init(ctx);
  
  res = siegel_modeq_eval(E, ctx, I, ell);
  if (res) modeq_all_isog_Q(nb_roots, all_M, nb_M, exp_array, E, ctx);
  
  modeq_clear(E);
  modeq_ctx_clear(ctx);
  return res;
}
