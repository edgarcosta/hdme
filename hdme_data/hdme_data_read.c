
#include "hdme_data.h"

void hdme_data_read(fmpq_mpoly_t pol, const char** vars, const char* name,
		    const fmpq_mpoly_ctx_t ctx)
{  
  char str[HDME_DATA_STR_LEN];
  char filename[HDME_DATA_FILE_LEN];
  char* success;
  FILE* file;
  int res;

  flint_sprintf(filename, "%s/%s", HDME_DATA_PATH, name);
  
  file = fopen(filename, "r");
  success = fgets(str, HDME_DATA_STR_LEN, file);
  if (success == NULL)
    {
      flint_printf("(mpoly_read) Error reading file %s\n", filename);
      fflush(stdout);
      flint_abort();
    }
  res = fmpq_mpoly_set_str_pretty(pol, str, vars, ctx);
  if (res == -1)
    {
      flint_printf("(mpoly_read) Error parsing string: %s\n", str);
      fflush(stdout);
      flint_abort();
    }
  fclose(file);
}
