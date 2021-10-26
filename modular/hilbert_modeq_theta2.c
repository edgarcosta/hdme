
#include "modular.h"

int hilbert_modeq_theta2(acb_ptr th2_vec, const acb_t t1, const acb_t t2,
			 const fmpz_poly_t beta, slong ell, slong delta, slong prec)
{
  acb_t z1, z2;
  acb_mat_t tau;
  fmpz_poly_mat_t m;
  sp2gz_t eta;
  arb_t tol;
  slong k;
  slong n = hilbert_nb_cosets(ell, delta);
  int res = 1;
  int v = MODEQ_VERBOSE;

  acb_init(z1);
  acb_init(z2);
  acb_mat_init(tau, 2, 2);
  fmpz_poly_mat_init(m, 2, 2);
  sp2gz_init(eta, 2);
  arb_init(tol);

  arb_one(tol);
  arb_mul_2exp_si(tol, tol, -MODEQ_RED_TOL_BITS);

  for (k = 0; k < n; k++)
    {
      if (v)
	{
	  flint_printf("(hilbert_modeq_theta2) Computing theta constants (%wd/%wd)\n", k+1, n);
	}
      if (res)
	{
	  hilbert_coset(m, k, ell, delta);
	  /* fmpz_poly_mat_print(m, "x"); */
	  hilbert_transform(z1, z2, m, t1, t2, delta, prec);
	  hilbert_scalar_div(z1, z2, beta, z1, z2, delta, prec);
	  hilbert_map(tau, z1, z2, delta, prec);
	  if (v)
	    {
	      flint_printf("(hilbert_modeq_theta2) Reduction... ");
	      fflush(stdout);
	    }
	  res = siegel_fundamental_domain(tau, eta, tau, tol, prec);
	  if (v) flint_printf("done.\n");
	}
      /* No star computation */
      if (res)
	{
	  res = theta2_unif(&th2_vec[16*k], tau, prec);
	}
      /* No renormalization necessary */
      if (!res)
	{
	  flint_printf("(siegel_modeq_theta) Warning: computation aborted due to low precision\n");
	  break;
	}
    }

  acb_clear(z1);
  acb_clear(z2);
  acb_mat_clear(tau);
  fmpz_poly_mat_clear(m);
  sp2gz_clear(eta);
  arb_clear(tol);
  return res;
}
