
#include "modular.h"

slong siegel_modeq_startprec_fmpq(fmpq* j, slong ell)
{
  slong h = siegel_modeq_height_fmpq(j);
  slong res = SIEGEL_START_PREC_MUL * n_pow(ell, 3) * n_clog(ell, 2)
    + (10 * siegel_nb_cosets(ell))/6 * h
    + SIEGEL_START_PREC_ADD;
  return 100 * (res/100 + 1);
}
