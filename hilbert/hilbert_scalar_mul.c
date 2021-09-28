
#include "hilbert.h"

void hilbert_scalar_mul(acb_t z1, acb_t z2, const fmpz_poly_t x, const acb_t t1,
			const acb_t t2, slong delta, slong prec)
{
  acb_t res1, res2;
  acb_t c;

  acb_init(res1);
  acb_init(res2);
  acb_init(c);

  hilbert_sigma1(c, x, delta, prec);
  acb_mul(res1, t1, c, prec);
  hilbert_sigma2(c, x, delta, prec);
  acb_mul(res2, t2, c, prec);

  acb_set(z1, res1);
  acb_set(z2, res2);
  
  acb_clear(res1);
  acb_clear(res2);
  acb_clear(c);
}
