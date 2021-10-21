
#include "modular.h"

slong hilbert_modeq_sym_igusa_startprec(fmpq* params, slong ell, slong len)
{
  slong h = hilbert_modeq_height(params, len);
  slong res = HILBERT_START_PREC_MUL * ell * (n_clog(ell, 2) + h);
  return 100 * (res/100 + 1);
}
