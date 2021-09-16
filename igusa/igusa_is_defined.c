
#include "igusa.h"

int igusa_is_defined(acb_srcptr j)
{
  return acb_is_finite(&j[0])
    && acb_is_finite(&j[1])
    && acb_is_finite(&j[2]);
}
