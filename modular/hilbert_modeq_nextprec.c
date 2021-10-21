
#include "modular.h"

slong hilbert_modeq_nextprec(slong current_prec)
{
  slong res = HILBERT_MUL_PREC * current_prec;
  return 100 * (res/100 + 1);
}
