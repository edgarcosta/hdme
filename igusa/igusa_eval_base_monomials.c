
#include "igusa.h"

void igusa_eval_base_monomials(acb_ptr ev, acb_srcptr I, slong wt, slong prec)
{
  fmpz_mpoly_ctx_t ctx;
  fmpz_mpoly_t mon;
  slong nb = igusa_nb_base_monomials(wt);
  slong k;

  cov_mpoly_ctx_init(ctx);
  fmpz_mpoly_init(mon, ctx);

  for (k = 0; k < nb; k++)
    {
      cov_base_monomial(mon, wt, k, ctx);
      cov_mpoly_eval(&ev[k], mon, I, ctx, prec);
    }

  fmpz_mpoly_clear(mon, ctx);
  cov_mpoly_ctx_clear(ctx);
}
