
#include "hilbert.h"

void hilbert_star(acb_t z, const fmpz_poly_mat_t m, const acb_t t1, const acb_t t2,
		  slong delta, slong prec)
{
  acb_t c, d, res, temp;
  acb_init(c);
  acb_init(d);
  acb_init(res);
  acb_init(temp);

  hilbert_sigma1(c, fmpz_poly_mat_entry(m, 1, 0), delta, prec);
  hilbert_sigma1(d, fmpz_poly_mat_entry(m, 1, 1), delta, prec);
  acb_mul(res, c, t1, prec);
  acb_add(res, res, d, prec);

  hilbert_sigma2(c, fmpz_poly_mat_entry(m, 1, 0), delta, prec);
  hilbert_sigma2(d, fmpz_poly_mat_entry(m, 1, 1), delta, prec);
  acb_mul(temp, c, t2, prec);
  acb_add(temp, temp, d, prec);
  acb_mul(res, res, temp, prec);

  acb_set(z, res);
  
  acb_clear(c);
  acb_clear(d);
  acb_clear(res);
  acb_clear(temp);
}
