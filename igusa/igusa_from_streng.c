
#include "igusa.h"

void igusa_from_streng(acb_ptr I, acb_srcptr S, slong prec)
{
  acb_mul_si(igusa_psi4(I), &S[0], 4);
  acb_mul_si(igusa_psi6(I), &S[1], 4);
  acb_mul_si(igusa_chi10(I), &S[2], -n_pow(2,12));
  acb_mul_si(igusa_chi12(I), &S[3], n_pow(2,15));
}
