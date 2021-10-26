
#include "modular.h"

void hilbert_modeq_igusa_C(acb_poly_t pol1, acb_poly_t pol2, acb_poly_t pol3,
			   acb_srcptr I_vec_beta, acb_srcptr I_vec_betabar,
			   slong ell, slong delta, slong prec)
{
  acb_ptr j1vec, j2vec, j3vec, j;
  acb_ptr xi, yi;
  slong n = hilbert_nb_cosets(ell, delta);
  slong k;
  int v = MODEQ_VERBOSE;

  j1vec = _acb_vec_init(2*n);
  j2vec = _acb_vec_init(2*n);
  j3vec = _acb_vec_init(2*n);
  xi = _acb_vec_init(2*n);
  yi = _acb_vec_init(2*n);
  j = _acb_vec_init(3);

  for (k = 0; k < n; k++)
    {
      igusa_from_cov(j, &I_vec_beta[4*k], prec);
      acb_set(&j1vec[k], &j[0]);
      acb_set(&j2vec[k], &j[1]);
      acb_set(&j3vec[k], &j[2]);
    }
  for (k = 0; k < n; k++)
    {
      igusa_from_cov(j, &I_vec_betabar[4*k], prec);
      acb_set(&j1vec[k+n], &j[0]);
      acb_set(&j2vec[k+n], &j[1]);
      acb_set(&j3vec[k+n], &j[2]);
    }

  for (k = 0; k < 2*n; k++)
    {
      acb_one(&xi[k]);
      acb_neg(&yi[k], &j1vec[k]);
    }

  if (v) flint_printf("(hilbert_modeq_sym_igusa_C) Building product trees...\n");
  product_tree_1(pol1, xi, yi, 2*n, prec);
  product_tree_2(pol2, xi, yi, j2vec, 2*n, prec);
  product_tree_2(pol3, xi, yi, j3vec, 2*n, prec);
  if (v) flint_printf("(hilbert_modeq_sym_igusa_C) Done.\n");

  _acb_vec_clear(j1vec, 2*n);
  _acb_vec_clear(j2vec, 2*n);
  _acb_vec_clear(j3vec, 2*n);
  _acb_vec_clear(xi, 2*n);
  _acb_vec_clear(yi, 2*n);
  _acb_vec_clear(j, 3);
}
