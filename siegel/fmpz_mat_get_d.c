
#include "siegel.h"

void fmpz_mat_get_d(fmpz_mat_t d, const fmpz_mat_t m)
{
  slong g = fmpz_mat_half_dim(m);
  slong j, k;
  for (j = 0; j < g; j++)
    {
      for (k = 0; k < g; k++)
        {
          fmpz_set(fmpz_mat_entry(d, j, k),
                   fmpz_mat_entry(m, j+g, k+g));
        }
    }
}
