
#include "modular.h"

int hilbert_modeq_igusa_eval_Q(fmpz_poly_struct* num_vec,
			       fmpz_t den, fmpq* rs, slong ell, slong delta)
{
  slong prec = hilbert_modeq_startprec(rs, ell, 2);
  acb_t r_acb, s_acb;
  acb_ptr I;
  acb_mat_t tau;
  acb_t t1, t2;
  sp2gz_t eta;
  fmpz_poly_t beta, betabar;
  acb_ptr th2_vec, I_vec_beta, I_vec_betabar;
  slong n = hilbert_nb_cosets(ell, delta);
  acb_poly_struct pol_vec_acb[3];
  int success = 0;
  int valid = 1;
  int v = MODEQ_VERBOSE;
  slong k;

  acb_init(r_acb);
  acb_init(s_acb);
  I = _acb_vec_init(4);
  acb_mat_init(tau, 2, 2);
  acb_init(t1);
  acb_init(t2);
  sp2gz_init(eta, 2);
  fmpz_poly_init(beta);
  fmpz_poly_init(betabar);
  th2_vec = _acb_vec_init(16*n);
  I_vec_beta = _acb_vec_init(4*n);
  I_vec_betabar = _acb_vec_init(4*n);
  for (k = 0; k < 3; k++) acb_poly_init(&pol_vec_acb[k]);

  valid = hilbert_splits(beta, ell, delta);
  if (valid) hilbert_conjugate(betabar, beta, delta);

  while (valid && !success && prec < MODEQ_MAX_PREC)
    {
      if (v) flint_printf("(hilbert_modeq_igusa_eval_Q) Start new run at precision %wd\n", prec);

      /* Do analytic computations */
      acb_set_fmpq(r_acb, &rs[0], prec);
      acb_set_fmpq(s_acb, &rs[1], prec);
      humbert_parametrize(I, r_acb, s_acb, delta, prec);
      if (v) flint_printf("(hilbert_modeq_igusa_eval_Q) Computing period matrix...\n");        
      valid = tau_from_igusa(tau, I, prec);
      if (v && valid) flint_printf("(hilbert_modeq_igusa_eval_Q) Done.\n");    
      if (v && !valid) flint_printf("(hilbert_modeq_igusa_eval_Q) Out of precision during computation of period matrix\n");
      if (valid)
	{
	  valid = hilbert_inverse(t1, t2, eta, tau, delta, prec);
	  if (v && !valid) flint_printf("(hilbert_modeq_igusa_eval_Q) Out of precision during inversion of Hilbert embedding\n");
	}
      if (valid) valid = hilbert_modeq_theta2(th2_vec, t1, t2, beta, ell, delta, prec);
      if (valid) modeq_cov(I_vec_beta, th2_vec, n, prec);
      if (valid) valid = hilbert_modeq_theta2(th2_vec, t1, t2, betabar, ell, delta, prec);
      if (valid) modeq_cov(I_vec_betabar, th2_vec, n, prec);
      
      if (v && !valid) flint_printf("(hilbert_modeq_igusa_eval_Q) Out of precision during complex computations\n");
      
      if (valid) /* Try rational reconstruction */
	{
	  hilbert_modeq_igusa_C(pol_vec_acb, I_vec_beta,
				I_vec_betabar, ell, delta, prec);
	  success = modeq_rational(num_vec, den, pol_vec_acb, 2*n, 3, prec);
	  if (v && !success) flint_printf("(hilbert_modeq_igusa_eval_Q) Not enough precision to recognize rational coefficients\n");
	  if (v && success) flint_printf("(hilbert_modeq_igusa_eval_Q) Heuristic success in recognizing coeffients: end of computation\n");
	}
      
      if (!success) /* Increase prec */
	{
	  prec = hilbert_modeq_nextprec(prec);
	  valid = 1;
	}
    }

  if (v && prec >= MODEQ_MAX_PREC) flint_printf("(hilbert_modeq_igusa_eval_Q) Maximum allowed precision reached: end of computation.\n");

  acb_clear(r_acb);
  acb_clear(s_acb);
  _acb_vec_clear(I, 4);
  acb_mat_clear(tau);
  acb_clear(t1);
  acb_clear(t2);
  sp2gz_clear(eta);
  fmpz_poly_clear(beta);
  fmpz_poly_clear(betabar);
  _acb_vec_clear(th2_vec, 16*n);
  _acb_vec_clear(I_vec_beta, 4*n);
  _acb_vec_clear(I_vec_betabar, 4*n);
  for (k = 0; k < 3; k++) acb_poly_clear(&pol_vec_acb[k]);

  return success;
}