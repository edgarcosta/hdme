
#include "modular.h"

int hilbert_modeq_eval(modeq_t R, modeq_ctx_t ctx, fmpz* I,
		       slong ell, slong delta)
{  
  hecke_t H;
  modeq_acb_t E;
  acb_t c;
  slong nb = 2*hilbert_nb_cosets(ell, delta);
  slong prec = hilbert_modeq_startprec(I, ell, delta);
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
      
      res = hecke_set_I_fmpz_hilbert(H, I, prec, delta);
      if (res) res = hecke_collect_hilbert_sym(H, ell, delta, prec);
      if (res) res = modeq_ctx_choose(ctx, hecke_all_I(H), nb, prec/div);
      if (res) modeq_product_trees(E, H, ctx, prec);
      if (res && (delta == 5))
	{
	  modeq_scalar(c, H, I, ctx, prec);
	  modeq_rescale(E, E, c, prec);
	  modeq_round(R, E);
	}
      else if (res) /* Delta is not 5 */
	{
	  res = modeq_rationalize(R, E, prec);
	}

      prec = siegel_modeq_nextprec(prec);
      stop = modeq_stop(res, prec);
    }
  
  hecke_clear(H);
  modeq_acb_clear(E);
  acb_clear(c);
  return res;
}
