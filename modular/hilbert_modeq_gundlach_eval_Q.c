
#include "modular.h"

int hilbert_modeq_gundlach_eval_Q(fmpz_poly_struct* num_vec,
				  fmpz_t den, fmpq* g, slong ell, slong delta)
{
  slong prec = hilbert_modeq_startprec(g, ell, 2);
  fmpz_poly_t beta, betabar;
  int stop = 0;
  int success = 1;
  int v = MODEQ_VERBOSE;
  slong k;
  slong n = hilbert_nb_cosets(ell, delta);
  int res;

  acb_ptr g_tau;
  acb_ptr j_tau;
  acb_ptr I_tau;
  acb_ptr th2_tau;
  acb_mat_t tau;
  fmpz_mat_t eta;
  acb_mat_t star;
  acb_t t1, t2;
  acb_ptr th2_vec;
  acb_ptr stardets;
  acb_ptr I_vec_beta;
  acb_ptr I_vec_betabar;
  acb_t scal;
  acb_poly_struct num_vec_acb[2];
  acb_t den_acb;
  fmpz_t rescale;
  acb_t rescale_acb;

  fmpz_poly_init(beta);
  fmpz_poly_init(betabar);
  g_tau = _acb_vec_init(2);
  j_tau = _acb_vec_init(3);
  I_tau = _acb_vec_init(4);
  th2_tau = _acb_vec_init(16);
  acb_mat_init(tau, 2, 2);
  fmpz_mat_init(eta, 4, 4);
  acb_mat_init(star, 2, 2);
  acb_init(t1);
  acb_init(t2);
  th2_vec = _acb_vec_init(2*16*n);
  stardets = _acb_vec_init(2*n);
  I_vec_beta = _acb_vec_init(4*n);
  I_vec_betabar = _acb_vec_init(4*n);
  acb_init(scal);
  for (k = 0; k < 2; k++) acb_poly_init(&num_vec_acb[k]);
  acb_init(den_acb);
  fmpz_init(rescale);
  acb_init(rescale_acb);
  
  success = hilbert_splits(beta, ell, delta);
  if (success) hilbert_conjugate(betabar, beta, delta);
  else
    {
      if (v) flint_printf("(hilbert_modeq_gundlach_eval_Q) Error: prime does not split correctly\n");
      stop = 1;
      res = 0;
    }
  
  hilbert_modeq_gundlach_rescale(rescale, g, ell, delta);
  acb_set_fmpz(rescale_acb, rescale);
  
  while (!stop)
    {
      if (v) flint_printf("(hilbert_modeq_gundlach_eval_Q) Start new run at precision %wd\n", prec);

      /* Do analytic computations */
      for (k = 0; k < 2; k++) acb_set_fmpq(&g_tau[k], &g[k], prec);
      igusa_from_gundlach(j_tau, g_tau, delta, prec);
      cov_from_igusa(I_tau, j_tau, prec);
      success = tau_theta2_from_igusa(tau, th2_tau, I_tau, prec);
      if (v && !success) flint_printf("(hilbert_modeq_gundlach_eval_Q) Out of precision when computing tau\n");
      if (success)
	{
	  success = hilbert_inverse(t1, t2, eta, tau, delta, prec);
	  if (v && !success) flint_printf("(hilbert_modeq_sym_igusa_eval_Q) Out of precision during inversion of Hilbert embedding\n");
	}      
      if (success)
	{
	  /* Rescale I_tau using eta */
	  cov_from_theta2(I_tau, th2_tau, prec);
	  siegel_star(star, eta, tau, prec);
	  acb_mat_det(scal, star, prec);
	  cov_rescale(I_tau, I_tau, scal, prec);
	  success = hilbert_modeq_theta2_star(th2_vec, stardets, t1, t2, beta,
					      ell, delta, prec)
	    && hilbert_modeq_theta2_star(&th2_vec[16*n], &stardets[n], t1, t2, betabar,
					 ell, delta, prec);	  
	  if (v && !success)
	    {
	      flint_printf("(hilbert_modeq_gundlach_eval_Q) Out of precision when computing theta constants\n");
	    }
	}
      if (success)
	{
	  modeq_cov(I_vec_beta, th2_vec, n, prec);
	  modeq_cov(I_vec_betabar, &th2_vec[16*n], n, prec);
	  hilbert_modeq_gundlach_scalar(scal, I_tau, stardets, ell, delta, prec);
	  hilbert_modeq_gundlach_num(num_vec_acb, I_vec_beta, I_vec_betabar,
				     scal, ell, delta, prec);
	  hilbert_modeq_gundlach_den(den_acb, I_vec_beta, I_vec_betabar,
				     scal, ell, delta, prec);
	  for (k = 0; k < 2; k++)
	    {
	      acb_poly_scalar_mul(&num_vec_acb[k], &num_vec_acb[k], rescale_acb, prec);
	    }	  
	  acb_mul_fmpz(den_acb, den_acb, rescale, prec);
	  success = modeq_round(num_vec, den, num_vec_acb, den_acb, 2*n, 2);
	  if (v && !success)
	    {
	      flint_printf("(hilbert_modeq_gundlach_eval_Q) Out of precision when recognizing integers\n");
	    }
	}      
      if (success)
	{
	  if (v) flint_printf("(hilbert_modeq_gundlach_eval_Q) Success at working precision %wd\n", prec);
	  modeq_simplify(num_vec, den, 2*n, 2);
	  stop = 1;
	  res = 1;
	}
      prec = hilbert_modeq_nextprec(prec);
      if (!stop && prec > MODEQ_MAX_PREC)
	{
	  flint_printf("(hilbert_modeq_gundlach_eval_Q) Reached maximal allowed precision %wd, abandon.\n", MODEQ_MAX_PREC);
	  stop = 1;
	  res = 0;
	}
    }
  
  fmpz_poly_clear(beta);
  fmpz_poly_clear(betabar);
  _acb_vec_clear(g_tau, 2);
  _acb_vec_clear(j_tau, 3);
  _acb_vec_clear(I_tau, 4);
  _acb_vec_clear(th2_tau, 16);
  acb_mat_clear(tau);
  fmpz_mat_clear(eta);
  acb_mat_clear(star);
  acb_clear(t1);
  acb_clear(t2);
  _acb_vec_clear(th2_vec, 2*16*n);
  _acb_vec_clear(stardets, 2*n);
  _acb_vec_clear(I_vec_beta, 4*n);
  _acb_vec_clear(I_vec_betabar, 4*n);
  acb_clear(scal);
  for (k = 0; k < 2; k++) acb_poly_clear(&num_vec_acb[k]);
  acb_clear(den_acb);
  fmpz_clear(rescale);
  acb_clear(rescale_acb);
  
  return res;      
}
