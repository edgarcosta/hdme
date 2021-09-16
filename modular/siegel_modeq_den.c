
#include "modular.h"

/* This uses the formula from the paper: "Evaluating modular equations
   in genus 2". We do not use etaR, since the c-block is zero and we
   know the determinant of the d-block. */

void siegel_modeq_den(acb_t den, acb_srcptr I_vec, const acb_t scal,
		      slong ell, slong prec)
{
  acb_t res, temp;
  slong k;

  acb_init(res);
  acb_init(temp);
  
  acb_one(res);
  for (k = 0; k < siegel_nb_cosets(ell); k++)
    {
      acb_sqr(temp, &I_vec[4*k+3], prec);
      acb_mul(res, res, temp, prec);
    }
  acb_mul(res, res, scal, prec);
  acb_set(den, res);

  acb_clear(res);
  acb_clear(temp);
}
