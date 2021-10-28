
#include "hilbert.h"

void humbert_parametrize(acb_ptr I, acb_srcptr rs, slong delta,
			 slong prec)
{
  acb_ptr AA1BB1B2;
  acb_t temp;

  AA1BB1B2 = _acb_vec_init(5);
  acb_init(temp);

  humbert_AA1BB1B2(AA1BB1B2, rs, delta, prec);
  if (delta == 33)
    {
      _acb_vec_set(I, AA1BB1B2, 4);
      acb_div(&I[0], &I[0], &AA1BB1B2[4], prec);
      acb_div(&I[2], &I[2], &AA1BB1B2[4], prec);
      acb_mul(temp, &I[0], &I[1], prec);
      acb_mul_si(&I[2], &I[2], -3, prec);
      acb_add(&I[2], &I[2], temp, prec);
      acb_div_si(&I[2], &I[2], 2, prec);
    }
  else
    {
      humbert_cov_from_AA1BB1B2(I, AA1BB1B2, prec);
    }

  _acb_vec_clear(AA1BB1B2, 5);
  acb_clear(temp);
}
