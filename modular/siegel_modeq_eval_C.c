
#include "modular.h"

int siegel_modeq_eval_C(acb_poly_struct* num_vec_0, acb_t den_0, acb_srcptr j_tau,
			slong ell, slong prec)
{
  slong n = siegel_nb_cosets(ell);

  acb_ptr I_tau;
  acb_mat_t tau;
  acb_ptr th2_vec;
  acb_ptr stardets;
  acb_ptr I_vec;
  acb_ptr th2_tau;
  acb_t scal;
  acb_t den;
  acb_poly_struct num_vec[3];
  slong k;
  int res = 1;
  int v = MODEQ_VERBOSE;

  I_tau = _acb_vec_init(4);
  acb_mat_init(tau, 2, 2);
  th2_vec = _acb_vec_init(16*n);
  stardets = _acb_vec_init(n);
  I_vec = _acb_vec_init(4*n);
  th2_tau = _acb_vec_init(16);
  acb_init(scal);
  acb_init(den);
  for (k = 0; k < 3; k++) acb_poly_init(&num_vec[k]);

  if (v) flint_printf("(siegel_modeq_eval_C) Run at precision %wd\n", prec);
  
  cov_from_igusa(I_tau, j_tau, prec);
  res = tau_theta2_from_igusa(tau, th2_tau, I_tau, prec);
  if (v && !res)
    {
      flint_printf("(siegel_modeq_eval_C) Out of precision when computing tau\n");
    }
  if (res)
    {
      res = siegel_modeq_theta2(th2_vec, stardets, tau, ell, prec);
      if (v && !res)
	{
	  flint_printf("(siegel_modeq_eval_C) Out of precision when computing theta constants at isogenous periods\n");
	}
    }
  if (res)
    {
      modeq_cov(I_vec, th2_vec, n, prec);
      cov_from_theta2(I_tau, th2_tau, prec);
      siegel_modeq_scalar(scal, I_tau, stardets, ell, prec);
      siegel_modeq_num(num_vec, I_vec, scal, ell, prec);
      siegel_modeq_den(den, I_vec, scal, ell, prec);
      if (v) flint_printf("(siegel_modeq_eval_C) Done.\n");
      for (k = 0; k < 3; k++) acb_poly_set(&num_vec_0[k], &num_vec[k]);
      acb_set(den_0, den);
    }
      
  _acb_vec_clear(I_tau, 4);
  acb_mat_clear(tau);
  _acb_vec_clear(th2_vec, 16*n);
  _acb_vec_clear(stardets, n);
  _acb_vec_clear(I_vec, 4*n);
  _acb_vec_clear(th2_tau, 16);
  acb_clear(scal);
  acb_clear(den);
  for (k = 0; k < 3; k++) acb_poly_clear(&num_vec[k]);
  return res;
}
