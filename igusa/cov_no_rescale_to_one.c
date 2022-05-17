
#include "igusa.h"

int cov_no_rescale_to_one(acb_srcptr I, slong prec)
{
  slong wt = 20;
  slong nb = cov_nb_base_monomials(wt);
  acb_ptr ev;
  slong i, j;
  int res = 1;

  ev = _acb_vec_init(nb);
  
  /* Compute weight 20 monomials and check overlap */
  cov_eval_base_monomials(ev, I, wt, prec);
  for (i = 0; i < nb; i++)
    {
      for (j = i+1; j < nb; j++)
	{
	  if (!acb_overlaps(&ev[i], &ev[j]))
	    {
	      res = 0;
	      break;
	    }
	}
      if (!res) break;
    }

  _acb_vec_clear(ev, nb);
  return res;
}
