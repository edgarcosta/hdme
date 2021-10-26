
#include "modular.h"

void hilbert_modeq_gundlach_num(acb_poly_struct* num_vec_acb,
				acb_srcptr I_vec_beta, acb_srcptr I_vec_betabar,
				const acb_t scal,
				slong ell, slong delta, slong prec)
{ 
  acb_ptr xi, yi, zi;
  acb_ptr G;
  acb_t temp;
  slong k;
  slong n = hilbert_nb_cosets(ell, delta);
  slong d = 2*n;
  
  xi = _acb_vec_init(d);
  yi = _acb_vec_init(d);
  zi = _acb_vec_init(d);
  acb_init(temp);
  G = _acb_vec_init(3);
  
  for (k = 0; k < n; k++)
    {
      gundlach_cov_from_igusa(G, &I_vec_beta[4*k], delta, prec);      
      acb_pow_si(&yi[k], &G[0], 5, prec);      
      acb_set(&xi[k], &G[2]);
      acb_sqr(&zi[k], &G[0], prec);
      acb_mul(&zi[k], &zi[k], &G[1], prec);
    }

  for (k = 0; k < n; k++)
    {
      gundlach_cov_from_igusa(G, &I_vec_betabar[4*k], delta, prec);      
      acb_pow_si(&yi[k+n], &G[0], 5, prec);      
      acb_set(&xi[k+n], &G[2]);
      acb_sqr(&zi[k+n], &G[0], prec);
      acb_mul(&zi[k+n], &zi[k+n], &G[1], prec);
    }
  
  flint_printf("(hilbert_modeq_gundlach_num) Building product trees...\n");
  product_tree_1(&num_vec_acb[0], xi, yi, d, prec);
  product_tree_2(&num_vec_acb[1], xi, yi, zi, d, prec);
  flint_printf("(hilbert_modeq_gundlach_num) Done.\n");

  for (k = 0; k < 2; k++)
    {
      acb_poly_scalar_mul(&num_vec_acb[k], &num_vec_acb[k], scal, prec);
    }

  _acb_vec_clear(xi, d);
  _acb_vec_clear(yi, d);
  _acb_vec_clear(zi, d);
  acb_clear(temp);
  _acb_vec_clear(G, 3);  
}
