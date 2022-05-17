
#include "siegel.h"

int fmpz_mat_is_scalar(const fmpz_mat_t m)
{
  slong n = fmpz_mat_nrows(m);
  slong j, k;
  
  if (n != fmpz_mat_ncols(m)) return 0;
  for (j = 0; j < n; j++)
    {
      for (k = 0; k < n; k++)
	{
	  if (j == k &&
	      !fmpz_equal(fmpz_mat_entry(m, j, k), fmpz_mat_entry(m, 0, 0)))
	    {
	      return 0;
	    }
	  if (j != k && !fmpz_is_zero(fmpz_mat_entry(m, j, k)))
	    {
	      return 0;
	    }
	}
    }
  return 1;  
}
