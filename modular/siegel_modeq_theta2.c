
#include "modular.h"

/* Storage: indices 16*k ... 16*(k+1)-1 contain theta constants
   computed at the image of tau under the k-th coset matrix, after
   reduction to the fundamental domain with some fixed
   tolerance. (Reduction matrix is stored as red[k]: we don't need
   it). */

int siegel_modeq_theta2(acb_ptr th2_vec, acb_ptr stardets,
			const acb_mat_t tau, slong ell, slong prec)
{
  slong k;
  acb_mat_t im;
  acb_mat_t red;
  acb_mat_t star;
  sp2gz_t eta;
  slong n = siegel_nb_cosets(ell);
  arb_t tol;
  int res;
  int v = SIEGEL_VERBOSE;

  acb_mat_init(im, 2, 2);
  acb_mat_init(red, 2, 2);
  acb_mat_init(star, 2, 2);
  sp2gz_init(eta, 2);
  arb_init(tol);

  arb_one(tol);
  arb_mul_2exp_si(tol, tol, -SIEGEL_RED_TOL_BITS);

  for (k = 0; k < n; k++)
    {
      if (v)
	{
	  flint_printf("(siegel_modeq_theta2) Computing theta constants (%wd/%wd)\n", k+1, n);
	}
      siegel_coset(eta, k, ell);
      res = siegel_transform(im, eta, tau, prec);
      if (res) res = siegel_fundamental_domain(red, eta, im, tol, prec);
      if (res)
	{
	  siegel_star(star, eta, im, prec);
	  acb_mat_det(&stardets[k], star, prec);
	  res = theta2_unif(&th2_vec[16*k], red, prec);
	}
      if (res) res = theta2_renormalize(&th2_vec[16*k], &th2_vec[16*k], prec);	
      if (!res)
	{
	  flint_printf("(siegel_modeq_theta) Warning: computation aborted due to low precision\n");
	  break;
	}
    }

  acb_mat_clear(im);
  acb_mat_clear(red);
  acb_mat_clear(star);
  sp2gz_clear(eta);
  arb_clear(tol);
  return res;
}
