
#include "igusa.h"

void cov_monomial(fmpz_mpoly_t mon, slong* exps, const fmpz_mpoly_ctx_t ctx)
{
  slong k;
  fmpz_mpoly_zero(mon);
  
  for (k = 0; k < 4; k++)
    {
      if (exps[k] < 0)
	{
	  flint_printf("(cov_monomial) Error: exponents must be nonnegative\n");
	  fflush(stdout);
	  flint_abort();
	}
    }
  fmpz_mpoly_set_coeff_ui_ui(mon, 1, (ulong*) exps, ctx);
}
