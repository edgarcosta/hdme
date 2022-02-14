
#include "modular.h"

int siegel_modeq_eval_Q(fmpz_poly_struct* num_vec,
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
  acb_t den_acb;
  acb_poly_struct num_vec_acb[3];
  fmpz_t rescale;
  acb_t rescale_acb;
  slong k;
  int stop = 0;
  int success;
  int res = 1;
  int v = MODEQ_VERBOSE;

  j_tau = _acb_vec_init(3);
  I_tau = _acb_vec_init(4);
  acb_mat_init(tau, 2, 2);
  th2_vec = _acb_vec_init(16*n);
  stardets = _acb_vec_init(n);
  I_vec = _acb_vec_init(4*n);
  th2_tau = _acb_vec_init(16);
  acb_init(scal);
  acb_init(den_acb);
  for (k = 0; k < 3; k++) acb_poly_init(&num_vec_acb[k]);
  fmpz_init(rescale);
  acb_init(rescale_acb);

  siegel_modeq_rescale(rescale, j, ell);
  acb_set_fmpz(rescale_acb, rescale);

  while (!stop)
    {
      if (v) flint_printf("(siegel_modeq_eval_Q) Start new run at precision %wd\n", prec);
      
      for (k = 0; k < 3; k++) acb_set_fmpq(&j_tau[k], &j[k], prec);
      cov_from_igusa(I_tau, j_tau, prec);
      success = tau_theta2_from_igusa(tau, th2_tau, I_tau, prec);
      if (v && !success)
	{
	  flint_printf("(siegel_modeq_eval_Q) Out of precision when computing tau\n");
	}
      if (success)
	{
	  success = siegel_modeq_theta2(th2_vec, stardets, tau, ell, prec);
	  if (v && !success)
	    {
	      flint_printf("(siegel_modeq_eval_Q) Out of precision when computing theta constants at isogenous periods\n");
	    }
	}
      if (success)
	{
	  modeq_cov(I_vec, th2_vec, n, prec);
	  cov_from_theta2(I_tau, th2_tau, prec);
	  siegel_modeq_scalar(scal, I_tau, stardets, ell, prec);
	  siegel_modeq_num(num_vec_acb, I_vec, scal, ell, prec);
	  siegel_modeq_den(den_acb, I_vec, scal, ell, prec);
	  for (k = 0; k < 3; k++)
	    {
	      acb_poly_scalar_mul(&num_vec_acb[k], &num_vec_acb[k],
				  rescale_acb, prec);
	    }
	  acb_mul_fmpz(den_acb, den_acb, rescale, prec);
	  success = modeq_round(num_vec, den,
				num_vec_acb, den_acb, n, 3);
	  if (v && !success)
	    {
	      flint_printf("(siegel_modeq_eval_Q) Out of precision when recognizing integers\n");
	      /* acb_printd(den_acb, 30); flint_printf("\n");
		 acb_printd(acb_poly_get_coeff_ptr(&num_vec_acb[2], 0), 30); flint_printf("\n"); */
	    }
	}
      if (success)
	{
	  if (v) flint_printf("(siegel_modeq_eval_Q) Success at working precision %wd\n", prec);
	  modeq_simplify(num_vec, den, n, 3);
	  stop = 1;
	}
      prec = siegel_modeq_nextprec(prec);
      if (!stop && prec > MODEQ_MAX_PREC)
	{
	  flint_printf("(siegel_modeq_eval_Q) Reached maximal allowed precision %wd, abandon.\n", MODEQ_MAX_PREC);
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
  acb_clear(den_acb);
  for (k = 0; k < 3; k++) acb_poly_clear(&num_vec_acb[k]);
  fmpz_clear(rescale);
  acb_clear(rescale_acb);
  return res;
}
