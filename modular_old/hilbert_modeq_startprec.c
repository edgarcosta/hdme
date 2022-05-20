
#include "modular.h"

slong hilbert_modeq_startprec(fmpq* params, slong ell, slong len)
{
  slong h = modeq_height_fmpq(params, len);
  slong res = HILBERT_START_PREC_MUL * ell * n_clog(ell, 2)
    + 1.2 * 10 * (ell+1)/3 * h
    + HILBERT_START_PREC_ADD;
  return 100 * (res/100 + 1);
}
