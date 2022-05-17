
#include "igusa.h"

void cov_monomial_exps(slong* exps, const fmpz_mpoly_t mon,
		       const fmpz_mpoly_ctx_t ctx)
{
  fmpz_mpoly_degrees_si(exps, mon, ctx);
}
