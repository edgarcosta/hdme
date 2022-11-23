
#include "hilbert.h"


void hilbert_get_mpoly(fmpq_mpoly_t pol, const char** vars, const char* name,
                       slong delta, const fmpq_mpoly_ctx_t ctx)
{
  char filename[HDME_DATA_FILE_LEN];
  flint_sprintf(filename, "hilbert/%wd/%s", delta, name);
  hdme_data_read(pol, vars, filename, ctx);
}
