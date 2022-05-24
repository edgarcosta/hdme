
#include "hecke.h"

void hecke_collect_print_status(int res, slong k, slong nb)
{
  if (!res)
    {
      flint_printf("(hecke_collect) Warning: computation aborted due to low precision\n");
    }
  else
    {
      if ((k+1) % 100 == 0)
	{
	  flint_printf("\n(hecke_collect) (%wd/%wd)", k+1, nb);
	}
      flint_printf("."); fflush(stdout);
    }
}
