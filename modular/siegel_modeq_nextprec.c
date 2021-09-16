
#include "modular.h"

slong siegel_modeq_nextprec(slong current_prec)
{
  slong res = SIEGEL_MUL_PREC * current_prec;
  return 100 * (res/100 + 1);
}
