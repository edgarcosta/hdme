
#include "hilbert.h"

void gundlach_get_mpoly(fmpq_mpoly_t pol, const char** vars, const char* name,
			slong delta, const fmpq_mpoly_ctx_t ctx)
{
  char str[HILBERT_MAX_STRLEN];
  char filename[HILBERT_MAX_STRLEN];
  char* success;
  FILE* file;
  int res;

  flint_sprintf(filename, "%s/%wd/%s", GUNDLACH_DATA_PATH, delta, name);
  /* flint_printf("(gundlach_get_mpoly) Reading %s\n", filename); */
  
  file = fopen(filename, "r");
  success = fgets(str, HILBERT_MAX_STRLEN, file);
  if (success == NULL)
    {
      flint_printf("(gundlach_get_mpoly) Error reading file %s\n", filename);
      fflush(stdout);
      flint_abort();
    }
  res = fmpq_mpoly_set_str_pretty(pol, str, vars, ctx);
  if (res == -1)
    {
      flint_printf("(gundlach_get_mpoly) Error parsing string: %s\n", str);
      fflush(stdout);
      flint_abort();
    }
  fclose(file);
}
