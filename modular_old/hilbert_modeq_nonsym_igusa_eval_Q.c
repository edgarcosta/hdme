
#include "modular.h"

int hilbert_modeq_nonsym_igusa_eval_Q(fmpz_poly_struct* num_vec,
				      fmpz_t den, fmpq* rs, slong ell, const fmpz_poly_t beta,
				      slong delta)
{
  slong prec = hilbert_modeq_startprec(rs, ell, 2);
  acb_ptr rs_acb;
  acb_ptr I;
  acb_mat_t tau;
  acb_ptr t;
  fmpz_mat_t eta;
  acb_ptr th2_vec, I_vec;
  slong n = hilbert_nb_cosets(ell, delta);
  acb_poly_struct pol_vec_acb[3];
  int success = 0;
  int valid = 1;
  int v = MODEQ_VERBOSE;
  slong k;

  rs_acb = _acb_vec_init(2);
  I = _acb_vec_init(4);
  acb_mat_init(tau, 2, 2);
  t = _acb_vec_init(2);
  fmpz_mat_init(eta, 4, 4);
  th2_vec = _acb_vec_init(16*n);
  I_vec = _acb_vec_init(4*n);
  for (k = 0; k < 3; k++)  acb_poly_init(&pol_vec_acb[k]);

  while (valid && !success && prec < MODEQ_MAX_PREC)
    {
      if (v) flint_printf("(hilbert_modeq_nonsym_igusa_eval_Q) Start new run at precision %wd\n", prec);

      /* Do analytic computations */
      acb_set_fmpq(&rs_acb[0], &rs[0], prec);
      acb_set_fmpq(&rs_acb[1], &rs[1], prec);
      hilbert_parametrize(I, rs_acb, delta, prec);
      if (v) flint_printf("(hilbert_modeq_nonsym_igusa_eval_Q) Computing period matrix...\n");        
      valid = tau_from_igusa(tau, I, prec);
      if (v && valid) flint_printf("(hilbert_modeq_nonsym_igusa_eval_Q) Done.\n");    
      if (v && !valid) flint_printf("(hilbert_modeq_nonsym_igusa_eval_Q) Out of precision during computation of period matrix\n");
      if (valid)
	{
	  valid = hilbert_inverse(t, eta, tau, delta, prec);
	  if (v && !valid) flint_printf("(hilbert_modeq_nonsym_igusa_eval_Q) Out of precision during inversion of Hilbert embedding\n");
	}
      if (valid) valid = hilbert_modeq_theta2(th2_vec, t, beta, ell, delta, prec);
      if (valid) modeq_cov(I_vec, th2_vec, n, prec);      
      if (v && !valid) flint_printf("(hilbert_modeq_nonsym_igusa_eval_Q) Out of precision during complex computations\n");
      
      if (valid) /* Try rational reconstruction */
	{
	  hilbert_modeq_nonsym_igusa_C(pol_vec_acb, I_vec, ell, delta, prec);
	  success = modeq_rational(num_vec, den, pol_vec_acb, n, 3, prec);
	  if (v && !success) flint_printf("(hilbert_modeq_nonsym_igusa_eval_Q) Not enough precision to recognize rational coefficients\n");
	  if (v && success) flint_printf("(hilbert_modeq_nonsym_igusa_eval_Q) Heuristic success in recognizing coeffients: end of computation\n");
	}
      
      if (!success)
	{
	  prec = hilbert_modeq_nextprec(prec);
	  valid = 1;
	}
    }

  if (v && prec >= MODEQ_MAX_PREC) flint_printf("(hilbert_modeq_nonsym_igusa_eval_Q) Maximum allowed precision reached: end of computation.\n");

  _acb_vec_clear(rs_acb, 2);
  _acb_vec_clear(I, 4);
  acb_mat_clear(tau);
  _acb_vec_clear(t, 2);
  fmpz_mat_clear(eta);
  _acb_vec_clear(th2_vec, 16*n);
  _acb_vec_clear(I_vec, 4*n);
  for (k = 0; k < 3; k++) acb_poly_clear(&pol_vec_acb[k]);

  return success;
}

