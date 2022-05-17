
#include "igusa.h"

void cov_mpoly_print(const fmpz_mpoly_t pol, const fmpz_mpoly_ctx_t ctx)
{
  char x[4][10] = {"I4", "I6prime", "I10", "I12"};  
  fmpz_mpoly_print_pretty(pol, x, ctx);
}
