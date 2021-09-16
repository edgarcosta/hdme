
#include "igusa.h"

void igusa_from_curve(acb_ptr j, const acb_poly_t crv, slong prec)
{
  acb_ptr I;
  
  I = _acb_vec_init(4);

  igusa_scalar_covariants(I, crv, prec);
  igusa_from_cov(j, I, prec);
  
  _acb_vec_clear(I, 4);
}
