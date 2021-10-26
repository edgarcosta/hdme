
#include "modular.h"

void hilbert_modeq_nonsym_igusa_C(acb_poly_struct* pol_vec, acb_srcptr I_vec, slong ell,
				  slong delta, slong prec)
{
  acb_ptr j1vec, j2vec, j3vec, j;
  acb_ptr xi, yi;
  slong n = hilbert_nb_cosets(ell, delta);
  slong k;
  int v = MODEQ_VERBOSE;

  j1vec = _acb_vec_init(n);
  j2vec = _acb_vec_init(n);
  j3vec = _acb_vec_init(n);
  xi = _acb_vec_init(n);
  yi = _acb_vec_init(n);
  j = _acb_vec_init(3);

  for (k = 0; k < n; k++)
    {
      igusa_from_cov(j, &I_vec[4*k], prec);
      acb_set(&j1vec[k], &j[0]);
      acb_set(&j2vec[k], &j[1]);
      acb_set(&j3vec[k], &j[2]);
    }

  for (k = 0; k < n; k++)
    {
      acb_one(&xi[k]);
      acb_neg(&yi[k], &j1vec[k]);
    }

  if (v) flint_printf("(hilbert_modeq_nonsym_igusa_C) Building product trees...\n");
  product_tree_1(&pol_vec[0], xi, yi, n, prec);
  product_tree_2(&pol_vec[1], xi, yi, j2vec, n, prec);
  product_tree_2(&pol_vec[2], xi, yi, j3vec, n, prec);
  if (v) flint_printf("(hilbert_modeq_nonsym_igusa_C) Done.\n");

  _acb_vec_clear(j1vec, n);
  _acb_vec_clear(j2vec, n);
  _acb_vec_clear(j3vec, n);
  _acb_vec_clear(xi, n);
  _acb_vec_clear(yi, n);
  _acb_vec_clear(j, 3);
}
