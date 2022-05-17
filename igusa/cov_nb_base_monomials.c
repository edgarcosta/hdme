
#include "igusa.h"

slong cov_nb_base_monomials(slong wt)
{
  if (w == 20) return 5;
  else if (w == 30) return 7;
  else if (w == 60) return 11;
  else
    {
      flint_printf("(cov_nb_base_monomials) Error: invalid weight %wd\n", wt);
      fflush(stdout);
      flint_abort();
    }
}
