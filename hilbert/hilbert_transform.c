
#include "hilbert.h"

void hilbert_transform(acb_t z1, acb_t z2, const fmpz_poly_mat_t m, const acb_t t1,
		       const acb_t t2, slong delta, slong prec)
{
  acb_t a, b, c, d, res1, res2, temp, sqrtd;
  acb_init(a);
  acb_init(b);
  acb_init(c);
  acb_init(d);
  acb_init(res1);
  acb_init(res2);
  acb_init(temp);
  acb_init(sqrtd);


  hilbert_sigma1(a, fmpz_poly_mat_entry(m, 0, 0), delta, prec);
  hilbert_sigma1(b, fmpz_poly_mat_entry(m, 0, 1), delta, prec);
  hilbert_sigma1(c, fmpz_poly_mat_entry(m, 1, 0), delta, prec);
  hilbert_sigma1(d, fmpz_poly_mat_entry(m, 1, 1), delta, prec);
  arb_sqrt_ui(acb_realref(sqrtd), delta, prec);
  acb_div(b, b, sqrtd, prec);
  acb_mul(c, c, sqrtd, prec);
  
  acb_mul(res1, a, t1, prec);
  acb_add(res1, res1, b, prec);
  acb_mul(temp, c, t1, prec);
  acb_add(temp, temp, d, prec);
  acb_div(res1, res1, temp, prec);

  hilbert_sigma2(a, fmpz_poly_mat_entry(m, 0, 0), delta, prec);
  hilbert_sigma2(b, fmpz_poly_mat_entry(m, 0, 1), delta, prec);
  hilbert_sigma2(c, fmpz_poly_mat_entry(m, 1, 0), delta, prec);
  hilbert_sigma2(d, fmpz_poly_mat_entry(m, 1, 1), delta, prec);
  acb_neg(sqrtd, sqrtd);
  acb_div(b, b, sqrtd, prec);
  acb_mul(c, c, sqrtd, prec);
  
  acb_mul(res2, a, t2, prec);
  acb_add(res2, res2, b, prec);
  acb_mul(temp, c, t2, prec);
  acb_add(temp, temp, d, prec);
  acb_div(res2, res2, temp, prec);

  acb_set(z1, res1);
  acb_set(z2, res2);

  acb_clear(a);
  acb_clear(b);
  acb_clear(c);
  acb_clear(d);
  acb_clear(res1);
  acb_clear(res2);
  acb_clear(temp);
  acb_clear(sqrtd);
}
