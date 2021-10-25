
#include "modular.h"

int hilbert_modeq_theta2_star(acb_ptr th2_vec, acb_ptr stardets,
			      const acb_t t1, const acb_t t2,
			      const fmpz_poly_t beta, slong ell, slong delta, slong prec)
{ 
  acb_t z1, z2;
  acb_mat_t im;
  acb_mat_t red;
  acb_mat_t star;
  acb_t hstar;
  fmpz_poly_mat_t m;
  sp2gz_t eta;
  arb_t tol;
  slong k;
  slong n = hilbert_nb_cosets(ell, delta);
  int res = 1;
  int v = HILBERT_VERBOSE;

  acb_init(z1);
  acb_init(z2);
  acb_mat_init(im, 2, 2);
  acb_mat_init(red, 2, 2);
  acb_mat_init(star, 2, 2);
  acb_init(hstar);
  fmpz_poly_mat_init(m, 2, 2);
  sp2gz_init(eta, 2);
  arb_init(tol);

  arb_one(tol);
  arb_mul_2exp_si(tol, tol, -SIEGEL_RED_TOL_BITS);

  for (k = 0; k < n; k++)
    {
      if (v)
	{
	  flint_printf("(hilbert_modeq_theta2_star) Computing theta constants (%wd/%wd)\n", k+1, n);
	}
      if (res)
	{
	  hilbert_coset(m, k, ell, delta);
	  /* fmpz_poly_mat_print(m, "x"); */
	  hilbert_transform(z1, z2, m, t1, t2, delta, prec);
	  hilbert_scalar_div(z1, z2, beta, z1, z2, delta, prec);
	  hilbert_map(im, z1, z2, delta, prec);
	  res = siegel_fundamental_domain(red, eta, im, tol, prec);
	}
      if (res)
	{
	  siegel_star(star, eta, im, prec);
	  acb_mat_det(&stardets[k], star, prec);
	  hilbert_star(hstar, m, t1, t2, delta, prec);
	  acb_mul(&stardets[k], &stardets[k], hstar, prec);
	  res = theta2_unif(&th2_vec[16*k], red, prec);
	}
      if (res) res = theta2_renormalize(&th2_vec[16*k], &th2_vec[16*k], prec);
      if (!res)
	{
	  flint_printf("(hilbert_modeq_theta2_star) Warning: computation aborted due to low precision\n");
	  break;
	}
    }

  acb_clear(z1);
  acb_clear(z2);
  acb_mat_clear(im);
  acb_mat_clear(red);
  acb_mat_clear(star);
  acb_clear(hstar);
  fmpz_poly_mat_clear(m);
  sp2gz_clear(eta);
  arb_clear(tol);
  return res;
}
