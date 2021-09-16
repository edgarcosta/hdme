
#include "modular.h"

void siegel_modeq_num(acb_poly_t num1, acb_poly_t num2, acb_poly_t num3,
		      acb_srcptr I_vec, const acb_t scal,
		      slong ell, slong prec)
{
  acb_ptr xi, yi, zi2, zi3;
  acb_t temp;
  slong k;
  slong d = siegel_nb_cosets(ell);
  
  xi = _acb_vec_init(d);
  yi = _acb_vec_init(d);
  zi2 = _acb_vec_init(d);
  zi3 = _acb_vec_init(d);
  acb_init(temp);

  for (k = 0; k < d; k++)
    {
      acb_mul(&yi[k], &I_vec[4*k+1], &I_vec[4*k+2], prec);
      acb_mul(&yi[k], &yi[k], &I_vec[4*k+3], prec);
      acb_neg(&yi[k], &yi[k]);
      
      acb_sqr(&xi[k], &I_vec[4*k+3], prec);

      acb_sqr(&zi2[k], &I_vec[4*k+1], prec);
      acb_mul(&zi2[k], &zi2[k], &I_vec[4*k], prec);
      acb_mul(&zi2[k], &zi2[k], &I_vec[4*k+3], prec);

      acb_pow_ui(&zi3[k], &I_vec[4*k+1], 5, prec);
    }
  flint_printf("(siegel_modeq_num) Building product trees...\n");
  product_tree_1(num1, xi, yi, d, prec);
  product_tree_2(num2, xi, yi, zi2, d, prec);
  product_tree_2(num3, xi, yi, zi3, d, prec);
  flint_printf("(siegel_modeq_num) Done.\n");
  
  acb_poly_scalar_mul(num1, num1, scal, prec);
  acb_poly_scalar_mul(num2, num2, scal, prec);
  acb_poly_scalar_mul(num3, num3, scal, prec);

  _acb_vec_clear(xi, d);
  _acb_vec_clear(yi, d);
  _acb_vec_clear(zi2, d);
  _acb_vec_clear(zi3, d);
  acb_clear(temp);
}
