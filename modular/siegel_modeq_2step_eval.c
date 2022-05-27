
#include "modular.h"

/* See also siegel_modeq_eval */

int siegel_modeq_2step_eval(modeq_t R, modeq_ctx_t ctx, fmpz* I, slong ell)
{  
  hecke_t H;
  modeq_acb_t E;
  acb_t c;
  slong nb = siegel_nb_T1_cosets(ell);
  slong prec = siegel_modeq_2step_startprec(I, ell);
  slong div = MODEQ_CTX_DIV_PREC;
  int res;
  int v = MODEQ_VERBOSE;
  int stop = 0;

  hecke_init(H, nb);
  modeq_acb_init(E);
  acb_init(c);
 
  while (!stop)
    {
      if (v) modeq_verbose_start(prec);
      
      res = hecke_set_I_fmpz(H, I, prec);
      if (res) res = hecke_collect_T1(H, ell, prec);
      if (res) res = modeq_ctx_choose(ctx, hecke_all_I(H), nb, prec/div);
      if (res)
	{
	  modeq_scalar(c, H, I, ctx, prec);
	  modeq_product_trees(E, H, ctx, prec);
	  modeq_rescale(E, E, c, prec);
	  res = modeq_round(R, E);
	}
      
      prec = siegel_modeq_nextprec(prec);
      stop = modeq_stop(res, prec);
    }
  
  hecke_clear(H);
  modeq_acb_clear(E);
  acb_clear(c);
  return res;
}
