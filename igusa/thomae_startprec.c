
#include "igusa.h"

slong thomae_startprec(slong prec)
{
  slong res = prec;
  while (res > 2 * THOMAE_LOWPREC - 2)
    {
      res = (res/2) + 1;
    }
  return res;
}
