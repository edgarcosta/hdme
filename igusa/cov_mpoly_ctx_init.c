
#include "igusa.h"

void cov_mpoly_ctx_init(fmpz_mpoly_ctx_t ctx)
{
  fmpz_mpoly_ctx_init(ctx, 4, ORD_LEX); 
}
