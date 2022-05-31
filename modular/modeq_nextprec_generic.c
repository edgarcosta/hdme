
#include "modular.h"

slong modeq_nextprec_generic(slong current_prec)
{
  slong res = MODEQ_MUL_PREC * current_prec;
  return 100 * (res/100 + 1);
}
