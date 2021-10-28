
#include "igusa.h"

int igusa_has_generic_automorphisms(acb_srcptr I, slong prec)
{
  acb_ptr ABCD;
  acb_t I6, R2;
  int res;

  ABCD = _acb_vec_init(4);
  acb_init(I6);
  acb_init(R2);

  igusa_I6(I6, I, prec);
  igusa_R2(R2, I, prec);
  igusa_clebsch(ABCD, I, prec);

  res = !acb_contains_zero(&I[3])
    && !acb_contains_zero(R2)
    && (!acb_contains_zero(&ABCD[0])
	|| !acb_contains_zero(&ABCD[1])
	|| !acb_contains_zero(&ABCD[2]));

  _acb_vec_clear(ABCD, 4);
  acb_clear(I6);
  acb_clear(R2);
  return res;
}
