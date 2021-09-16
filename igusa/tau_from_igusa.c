
#include "igusa.h"

int tau_from_igusa(acb_mat_t tau, acb_srcptr I, slong prec)
{
  slong perm, signs;
  acb_poly_t crv;
  acb_ptr roots;
  acb_ptr ros;
  acb_ptr th2, th4;
  int res;

  acb_poly_init(crv);
  roots = _acb_vec_init(6);
  ros = _acb_vec_init(3);
  th2 = _acb_vec_init(16);
  th4 = _acb_vec_init(16);

  res = mestre(crv, I, prec);
  if (res) res = thomae_roots(roots, crv, prec);
  if (res) res = thomae_correct_signs(&perm, &signs, roots, I, prec);
  if (res)
    {
      flint_printf("(tau_from_igusa) Computing period matrix at high precision...\n");
      thomae_reorder(roots, roots, perm);
      thomae_rosenhain(ros, roots, prec);
      thomae_theta4(th4, ros, prec);
      thomae_theta2(th2, th4, ros, signs, prec);
      res = theta2_inverse(tau, th2, prec);
      flint_printf("(tau_from_igusa) Done.\n");
    }

  acb_poly_clear(crv);
  _acb_vec_clear(roots, 6);
  _acb_vec_clear(ros, 3);
  _acb_vec_clear(th2, 16);
  _acb_vec_clear(th4, 16);
  return res;
}
