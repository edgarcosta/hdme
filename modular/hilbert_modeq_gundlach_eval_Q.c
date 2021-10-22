
#include "modular.h"

int hilbert_modeq_gundlach_eval_Q(fmpz_poly_t num1, fmpz_poly_t num2,
				  fmpz_t den, fmpq* g, slong ell, slong delta)
{
  slong prec = hilbert_modeq_sym_igusa_startprec(g, ell, 2);
  fmpz_poly_t beta, betabar;
  int stop = 0;
  int success = 1;
  int v = HILBERT_VERBOSE;
  slong k;
  slong n = hilbert_nb_cosets(ell);
  int res;

  acb_ptr g_tau;
  acb_ptr j_tau;
  acb_ptr I_tau;
  acb_ptr th2_tau;
  sp2gz_t eta;
  acb_t t1, t2;
  acb_ptr th2_vec;
  acb_ptr stardets;
  acb_ptr I_vec_beta;
  acb_ptr I_vec_betabar;
  acb_t scal;
  acb_poly_t num1_acb, num2_acb;
  acb_t den_acb;
  fmpz_t rescale;
  acb_t rescale_acb;

  fmpz_poly_init(beta);
  fmpz_poly_init(betabar);
  g_tau = _acb_vec_init(2);
  j_tau = _acb_vec_init(3);
  I_tau = _acb_vec_init(4);
  th2_tau = _acb_vec_init(16);
  sp2gz_init(eta, 2);
  acb_init(t1);
  acb_init(t2);
  th2_vec = _acb_vec_init(2*16*n);
  stardets = _acb_vec_init(2*n);
  I_vec_beta = _acb_vec_init(4*n);
  I_vec_betabar = _acb_vec_init(4*n);
  acb_init(scal);
  acb_poly_init(num1_acb);
  acb_poly_init(num2_acb);
  acb_init(den_acb);
  fmpz_init(rescale);
  acb_init(rescale_acb);
  
  valid = hilbert_splits(beta, ell, delta);
  if (valid) hilbert_conjugate(betabar, beta, delta);
  
  hilbert_modeq_gundlach_fmpq_rescale(rescale, g, ell, delta);
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
	  valid = hilbert_inverse(t1, t2, eta, tau, delta, prec);
	  if (v && !success) flint_printf("(hilbert_modeq_sym_igusa_eval_Q) Out of precision during inversion of Hilbert embedding\n");
	}      
      if (success)
	{
	  success = hilbert_modeq_theta2_star(th2_vec, stardets, t1, t2, beta,
					      ell, delta, prec)
	    && hilbert_modeq_theta2_star(&th2vec[16*n], &stardets[n], t1, t2, betabar,
					 ell, delta, prec);	  
	  if (!success)
	    {
	      flint_printf("(hilbert_modeq_gundlach_eval_Q) Out of precision when computing theta constants\n");
	    }
	}
      /* Need to rescale I_tau using eta. */
      if (success)
	{
	  hilbert_modeq_cov(I_vec_beta, th2_vec, ell, delta, prec);
	  hilbert_modeq_cov(I_vec_betabar, &th2_vec[16*n], ell, delta, prec);
	  hilbert_modeq_gundlach_scalar(scal, I_tau, stardets, ell, delta, prec);
	  hilbert_modeq_gundlach_num(num1, num2, I_vec_beta, I_vec_betabar,
				     scal, ell, delta, prec);
	  hilbert_modeq_gundlach_den(den, I_vec_beta, I_vec_betabar,
				     scal, ell, delta, prec);
	  acb_poly_scalar_mul(num1_acb, num1_acb, rescale_acb, prec);
	  acb_poly_scalar_mul(num2_acb, num2_acb, rescale_acb, prec);
	  acb_mul_fmpz(den_acb, den_acb, rescale, prec);
	  success = hilbert_modeq_gundlach_round(num1, num2, den,
						 num1_acb, num2_acb, den_acb,
						 ell, delta);
	  if (!success)
	    {
	      flint_printf("(hilbert_modeq_gundlach_eval_Q) Out of precision when recognizing integers\n");
	    }
	}      
      if (success)
	{
	  flint_printf("(hilbert_modeq_gundlach_eval_Q) Success at working precision %wd\n", prec);
	  hilbert_modeq_gundlach_simplify(num1, num2, den, ell, delta);
	  stop = 1;
	}
      prec = hilbert_modeq_nextprec(prec);
      if (!stop && prec > HILBERT_MAX_PREC)
	{
	  flint_printf("(hilbert_modeq_gundlach_eval_Q) Reached maximal allowed precision %wd, abandon.\n", HILBERT_MAX_PREC);
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
  sp2gz_clear(eta);
  acb_clear(t1);
  acb_clear(t2);
  _acb_vec_clear(th2_vec, 2*16*n);
  _acb_vec_clear(stardets, 2*n);
  _acb_vec_clear(I_vec_beta, 4*n);
  _acb_vec_clear(I_vec_betabar, 4*n);
  acb_clear(scal);
  acb_poly_clear(num1_acb);
  acb_poly_clear(num2_acb);
  acb_clear(den_acb);
  fmpz_clear(rescale);
  acb_clear(rescale_acb);
  
  return res;      
}
