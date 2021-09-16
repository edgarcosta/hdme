
#include "modular.h"

int siegel_modeq_eval_Q(fmpz_poly_t num1, fmpz_poly_t num2, fmpz_poly_t num3,
			 fmpz_t den, fmpq* j, slong ell)
{
  slong prec = siegel_modeq_startprec_fmpq(j, ell);
  slong n = siegel_nb_cosets(ell);

  acb_ptr j_tau;
  acb_ptr I_tau;
  acb_mat_t tau;
  acb_ptr th2_vec;
  acb_ptr stardets;
  acb_ptr I_vec;
  acb_ptr th2_tau;
  acb_t scal;
  acb_poly_t num1_acb, num2_acb, num3_acb;
  acb_t den_acb;
  fmpz_t rescale;
  acb_t rescale_acb;
  slong k;
  int stop = 0;
  int success;
  int res = 1;

  j_tau = _acb_vec_init(3);
  I_tau = _acb_vec_init(4);
  acb_mat_init(tau, 2, 2);
  th2_vec = _acb_vec_init(16*n);
  stardets = _acb_vec_init(n);
  I_vec = _acb_vec_init(4*n);
  th2_tau = _acb_vec_init(16);
  acb_init(scal);
  acb_poly_init(num1_acb);
  acb_poly_init(num2_acb);
  acb_poly_init(num3_acb);
  acb_init(den_acb);
  fmpz_init(rescale);
  acb_init(rescale_acb);

  siegel_modeq_fmpq_rescale(rescale, j, ell);
  acb_set_fmpz(rescale_acb, rescale);

  while (!stop)
    {
      flint_printf("(siegel_modeq_eval_Q) Start new run at precision %wd\n", prec);
      
      for (k = 0; k < 3; k++) acb_set_fmpq(&j_tau[k], &j[k], prec);
      cov_from_igusa(I_tau, j_tau, prec);
      success = tau_theta2_from_igusa(tau, th2_tau, I_tau, prec);
      if (!success)
	{
	  flint_printf("(siegel_modeq_eval_Q) Out of precision when computing tau\n");
	}
      if (success)
	{
	  success = siegel_modeq_theta2(th2_vec, stardets, tau, ell, prec);
	  if (!success)
	    {
	      flint_printf("(siegel_modeq_eval_Q) Out of precision when computing theta constants at isogenous periods\n");
	    }
	}
      if (success)
	{
	  siegel_modeq_cov(I_vec, th2_vec, ell, prec);
	  cov_from_theta2(I_tau, th2_tau, prec);
	  siegel_modeq_scalar(scal, I_tau, stardets, ell, prec);
	  siegel_modeq_num(num1_acb, num2_acb, num3_acb, I_vec, scal, ell, prec);
	  siegel_modeq_den(den_acb, I_vec, scal, ell, prec);
	  acb_poly_scalar_mul(num1_acb, num1_acb, rescale_acb, prec);
	  acb_poly_scalar_mul(num2_acb, num2_acb, rescale_acb, prec);
	  acb_poly_scalar_mul(num3_acb, num3_acb, rescale_acb, prec);
	  acb_mul_fmpz(den_acb, den_acb, rescale, prec);
	  success = siegel_modeq_round(num1, num2, num3, den,
				       num1_acb, num2_acb, num3_acb, den_acb, ell);
	  if (!success)
	    {
	      flint_printf("(siegel_modeq_eval_Q) Out of precision when recognizing integers\n");
	    }
	}
      if (success)
	{
	  flint_printf("(siegel_modeq_eval_Q) Success at working precision %wd\n", prec);
	  siegel_modeq_simplify(num1, num2, num3, den, ell);
	  stop = 1;
	}
      prec = siegel_modeq_nextprec(prec);
      if (!stop && prec > SIEGEL_MAX_PREC)
	{
	  flint_printf("(siegel_modeq_eval_Q) Reached maximal allowed precision %wd, abandon.\n", SIEGEL_MAX_PREC);
	  stop = 1;
	  res = 0;
	}
    }
      
  _acb_vec_clear(j_tau, 3);
  _acb_vec_clear(I_tau, 4);
  acb_mat_clear(tau);
  _acb_vec_clear(th2_vec, 16*n);
  _acb_vec_clear(stardets, n);
  _acb_vec_clear(I_vec, 4*n);
  _acb_vec_clear(th2_tau, 16);
  acb_clear(scal);
  acb_poly_clear(num1_acb);
  acb_poly_clear(num2_acb);
  acb_poly_clear(num3_acb);
  acb_clear(den_acb);
  fmpz_clear(rescale);
  acb_clear(rescale_acb);
  return res;
}
