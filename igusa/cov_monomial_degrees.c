
#include "igusa.h"

void cov_monomial_degrees(slong* exps, const fmpz_mpoly_t mon,
			  const fmpz_mpoly_ctx_t ctx)
{
  fmpz_mpoly_degree_si(exps, mon, ctx);
}
