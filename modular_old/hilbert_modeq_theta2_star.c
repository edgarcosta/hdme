
#include "modular.h"

int hilbert_modeq_theta2_star(acb_ptr th2_vec, acb_ptr stardets,
			      acb_srcptr t,
			      const fmpz_poly_t beta, slong ell, slong delta, slong prec)
{
  acb_ptr z;
  acb_mat_t im;
  acb_mat_t red;
  acb_mat_t star;
  acb_t hstar;
  fmpz_poly_mat_t m;
  fmpz_mat_t eta;
  arb_t tol;
  slong k;
  slong n = hilbert_nb_cosets(ell, delta);
  int res = 1;
  int v = MODEQ_VERBOSE;

  z = _acb_vec_init(2);
  acb_mat_init(im, 2, 2);
  acb_mat_init(red, 2, 2);
  acb_mat_init(star, 2, 2);
  acb_init(hstar);
  fmpz_poly_mat_init(m, 2, 2);
  fmpz_mat_init(eta, 4, 4);
  arb_init(tol);

  arb_one(tol);
  arb_mul_2exp_si(tol, tol, -MODEQ_RED_TOL_BITS);

  if(v) flint_printf("(hilbert_modeq_theta2_star) Computing theta constants (%wd)", n);
  fflush(stdout);
  for (k = 0; k < n; k++)
    {
      if (v)
	{	  
	  if ((k+1) % 100 == 0)
	    {
	      flint_printf("\n(hilbert_modeq_theta2_star) (%wd/%wd)", k+1, n);
	    }
	  flint_printf("."); fflush(stdout);
	}
      if (res)
	{
	  hilbert_coset(m, k, ell, delta);
	  /* fmpz_poly_mat_print(m, "x"); */
	  hilbert_transform(z, m, t, delta, prec);
	  hilbert_scalar_div(z, beta, z, delta, prec);
	  hilbert_map(im, z, delta, prec);
	  if (v)
	    {
	      /*flint_printf("(hilbert_modeq_theta2_star) Reduction... ");
		fflush(stdout);*/
	    }
	  res = siegel_fundamental_domain(red, eta, im, tol, prec);
	  /* if (v) flint_printf("done.\n"); */
	}
	
      if (res)
	{
	  siegel_star(star, eta, im, prec);
	  acb_mat_det(&stardets[k], star, prec);
	  hilbert_star(hstar, m, t, delta, prec);
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
  if (v) flint_printf("\n");

  _acb_vec_clear(z, 2);
  acb_mat_clear(im);
  acb_mat_clear(red);
  acb_mat_clear(star);
  acb_clear(hstar);
  fmpz_poly_mat_clear(m);
  fmpz_mat_clear(eta);
  arb_clear(tol);
  return res;
}
