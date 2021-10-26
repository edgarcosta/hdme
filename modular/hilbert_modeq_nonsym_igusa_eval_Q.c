
#include "modular.h"

int hilbert_modeq_nonsym_igusa_eval_Q(fmpz_poly_t num1, fmpz_poly_t num2, fmpz_poly_t num3,
				      fmpz_t den, fmpq* rs, slong ell, const fmpz_poly_t beta,
				      slong delta)
{
  slong prec = hilbert_modeq_startprec(rs, ell, 2);
  acb_t r_acb, s_acb;
  acb_ptr I;
  acb_mat_t tau;
  acb_t t1, t2;
  sp2gz_t eta;
  acb_ptr th2_vec, I_vec;
  fmpq_poly_t pol1, pol2, pol3;
  fmpz_t den1, den2, den3;
  slong n = hilbert_nb_cosets(ell, delta);
  acb_poly_t pol1_acb, pol2_acb, pol3_acb;
  int success = 0;
  int valid = 1;
  int v = MODEQ_VERBOSE;

  acb_init(r_acb);
  acb_init(s_acb);
  I = _acb_vec_init(4);
  acb_mat_init(tau, 2, 2);
  acb_init(t1);
  acb_init(t2);
  sp2gz_init(eta, 2);
  th2_vec = _acb_vec_init(16*n);
  I_vec = _acb_vec_init(4*n);
  fmpq_poly_init(pol1);
  fmpq_poly_init(pol2);
  fmpq_poly_init(pol3);
  fmpz_init(den1);
  fmpz_init(den2);
  fmpz_init(den3);
  acb_poly_init(pol1_acb);
  acb_poly_init(pol2_acb);
  acb_poly_init(pol3_acb);

  while (valid && !success && prec < MODEQ_MAX_PREC)
    {
      if (v) flint_printf("(hilbert_modeq_nonsym_igusa_eval_Q) Start new run at precision %wd\n", prec);

      /* Do analytic computations */
      acb_set_fmpq(r_acb, &rs[0], prec);
      acb_set_fmpq(s_acb, &rs[1], prec);
      hilbert_parametrize(I, r_acb, s_acb, delta, prec);
      if (v) flint_printf("(hilbert_modeq_nonsym_igusa_eval_Q) Computing period matrix...\n");        
      valid = tau_from_igusa(tau, I, prec);
      if (v && valid) flint_printf("(hilbert_modeq_nonsym_igusa_eval_Q) Done.\n");    
      if (v && !valid) flint_printf("(hilbert_modeq_nonsym_igusa_eval_Q) Out of precision during computation of period matrix\n");
      if (valid)
	{
	  valid = hilbert_inverse(t1, t2, eta, tau, delta, prec);
	  if (v && !valid) flint_printf("(hilbert_modeq_nonsym_igusa_eval_Q) Out of precision during inversion of Hilbert embedding\n");
	}
      if (valid) valid = hilbert_modeq_theta2(th2_vec, t1, t2, beta, ell, delta, prec);
      if (valid) modeq_cov(I_vec, th2_vec, n, prec);      
      if (v && !valid) flint_printf("(hilbert_modeq_nonsym_igusa_eval_Q) Out of precision during complex computations\n");
      
      if (valid) /* Try rational reconstruction */
	{
	  hilbert_modeq_nonsym_igusa_C(pol1_acb, pol2_acb, pol3_acb, I_vec,
				       ell, delta, prec);
	  success = modeq_rational_poly(pol1, pol1_acb, n, prec);
	  if (success) success = modeq_rational_poly(pol2, pol2_acb, n-1, prec);
	  if (success) success = modeq_rational_poly(pol3, pol3_acb, n-1, prec);
	  if (v && !success) flint_printf("(hilbert_modeq_nonsym_igusa_eval_Q) Not enough precision to recognize rational coefficients\n");
	  if (v && success) flint_printf("(hilbert_modeq_nonsym_igusa_eval_Q) Heuristic success in recognizing coeffients: end of computation\n");
	}
      
      if (success) /* Set output values */
	{
	  fmpq_poly_get_denominator(den1, pol1);
	  fmpq_poly_get_denominator(den2, pol2);
	  fmpq_poly_get_denominator(den3, pol3);
	  fmpz_lcm(den, den1, den2);
	  fmpz_lcm(den, den, den3);
	  fmpq_poly_scalar_mul_fmpz(pol1, pol1, den);
	  fmpq_poly_scalar_mul_fmpz(pol2, pol2, den);
	  fmpq_poly_scalar_mul_fmpz(pol3, pol3, den);
	  fmpq_poly_get_numerator(num1, pol1);
	  fmpq_poly_get_numerator(num2, pol2);
	  fmpq_poly_get_numerator(num3, pol3);
	}
      else
	{
	  prec = hilbert_modeq_nextprec(prec);
	  valid = 1;
	}
    }

  if (v && prec >= MODEQ_MAX_PREC) flint_printf("(hilbert_modeq_nonsym_igusa_eval_Q) Maximum allowed precision reached: end of computation.\n");

  acb_clear(r_acb);
  acb_clear(s_acb);
  _acb_vec_clear(I, 4);
  acb_mat_clear(tau);
  acb_clear(t1);
  acb_clear(t2);
  sp2gz_clear(eta);
  _acb_vec_clear(th2_vec, 16*n);
  _acb_vec_clear(I_vec, 4*n);
  fmpq_poly_clear(pol1);
  fmpq_poly_clear(pol2);
  fmpq_poly_clear(pol3);
  fmpz_clear(den1);
  fmpz_clear(den2);
  fmpz_clear(den3);
  acb_poly_clear(pol1_acb);
  acb_poly_clear(pol2_acb);
  acb_poly_clear(pol3_acb);

  return success;
}

