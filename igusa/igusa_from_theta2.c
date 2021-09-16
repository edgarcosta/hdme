
#include "igusa.h"

void igusa_from_theta2(acb_ptr j, acb_srcptr theta2, slong prec)
{
  acb_ptr h;
  acb_t h10inv;
    
  h = _acb_vec_init(4);
  acb_init(h10inv);

  igusa_h(h, theta2, prec);
  
  acb_inv(h10inv, &h[2], prec);
  acb_mul(&j[0], &h[0], &h[1], prec);
  acb_mul(&j[0], &j[0], h10inv, prec);

  acb_sqr(h10inv, h10inv, prec); /* Now 1/h10^2 */
  acb_sqr(&j[1], &h[0], prec);
  acb_mul(&j[1], &j[1], &h[3], prec);
  acb_mul(&j[1], &j[1], h10inv, prec);
  
  acb_pow_ui(&j[2], &h[0], 5, prec);
  acb_mul(&j[2], &j[2], h10inv, prec);

  _acb_vec_clear(h, 4);
  acb_clear(h10inv);
}
