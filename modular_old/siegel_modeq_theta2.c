
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
  fmpz_mat_t eta;
  slong n = siegel_nb_cosets(ell);
  arb_t tol;
  int res;
  int v = MODEQ_VERBOSE;
  /* slong i;*/

  acb_mat_init(im, 2, 2);
  acb_mat_init(red, 2, 2);
  acb_mat_init(star, 2, 2);
  fmpz_mat_init(eta, 4, 4);
  arb_init(tol);

  arb_one(tol);
  arb_mul_2exp_si(tol, tol, -MODEQ_RED_TOL_BITS);

  if (v) flint_printf("(siegel_modeq_theta2) Computing theta constants (%wd)", n);
  fflush(stdout);
  for (k = 0; k < n; k++)
    {
      if (v)
	{
	  if ((k+1) % 100 == 0)
	    {
	      flint_printf("\n(siegel_modeq_theta2) (%wd/%wd)", k+1, n);
	    }
	  flint_printf("."); fflush(stdout);
	}
      siegel_coset(eta, k, ell);
      /*fmpz_mat_print_pretty(eta);*/
      res = siegel_transform(im, eta, tau, prec);
      if (res)
	{
	  res = siegel_fundamental_domain(red, eta, im, tol, prec);
	  /* acb_mat_printd(tau, 10);
	     acb_mat_printd(im, 10);
	     acb_mat_printd(red, 10); */
	}
      if (res)
	{
	  siegel_star(star, eta, im, prec);
	  acb_mat_det(&stardets[k], star, prec);
	  res = theta2_unif(&th2_vec[16*k], red, prec);
	}
      if (res)
	{
	  res = theta2_renormalize(&th2_vec[16*k], &th2_vec[16*k], prec);
	  /* for (i = 0; i < 16; i++)
	    {
	      acb_printd(&th2_vec[16*k+i], 10); flint_printf("\n");
	      } */
	}
      if (!res)
	{
	  flint_printf("(siegel_modeq_theta2) Warning: computation aborted due to low precision\n");
	  break;
	}
    }
  if (v) flint_printf("\n");

  acb_mat_clear(im);
  acb_mat_clear(red);
  acb_mat_clear(star);
  fmpz_mat_clear(eta);
  arb_clear(tol);
  return res;
}
