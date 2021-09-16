
#include "acb_mat_extras.h"

void
acb_mat_set_arb_arb(acb_mat_t z, const arb_mat_t re, const arb_mat_t im)
{
  slong nrows = acb_mat_nrows(re);
  slong ncols = acb_mat_ncols(re);
  slong i, j;

  for (i = 0; i < nrows; i++)
    {
      for (j = 0; j < ncols; j++)
	{
	  acb_set_arb_arb(acb_mat_entry(z, i, j),
			  arb_mat_entry(re, i, j), arb_mat_entry(im, i, j));
	}
    }
}
