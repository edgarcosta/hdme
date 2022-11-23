
#include "theta.h"

int theta_0123_der_naive(acb_ptr th, acb_mat_t dth,
			 const acb_mat_t tau, slong prec)
{
  acb_mat_t tau_pert;
  arb_t eps;
  acb_ptr value, value_pert;
  mag_t error;
  slong k;
  slong n = 4;
  int res = 1;

  acb_mat_init(tau_pert, 2, 2);
  arb_init(eps);
  mag_init(error);
  value = _acb_vec_init(n);
  value_pert = _acb_vec_init(n);

  theta_der_set_pert(eps, prec);
  res = theta_der_set_error(error, tau, prec);
  theta_0123_naive(value, tau, prec);

  acb_mat_set(tau_pert, tau);
  acb_add_arb(acb_mat_entry(tau_pert, 0, 0),
	      acb_mat_entry(tau_pert, 0, 0), eps, prec);
  res = res && theta_0123_naive(value_pert, tau_pert, prec);
  _acb_vec_sub(value_pert, value_pert, value, n, prec);
  _acb_vec_scalar_div_arb(value_pert, value_pert, n, eps, prec);
  for (k = 0; k < n; k++)
    {
      acb_set(acb_mat_entry(dth, k, 0), &value_pert[k]);
    }

  acb_mat_set(tau_pert, tau);
  acb_add_arb(acb_mat_entry(tau_pert, 1, 1),
	      acb_mat_entry(tau_pert, 1, 1), eps, prec);
  res = res && theta_0123_naive(value_pert, tau_pert, prec);
  _acb_vec_sub(value_pert, value_pert, value, n, prec);
  _acb_vec_scalar_div_arb(value_pert, value_pert, n, eps, prec);
  for (k = 0; k < n; k++)
    {
      acb_set(acb_mat_entry(dth, k, 1), &value_pert[k]);
    }

  acb_mat_set(tau_pert, tau);
  acb_add_arb(acb_mat_entry(tau_pert, 0, 1),
	      acb_mat_entry(tau_pert, 0, 1), eps, prec);
  acb_add_arb(acb_mat_entry(tau_pert, 1, 0),
	      acb_mat_entry(tau_pert, 1, 0), eps, prec);
  res = res && theta_0123_naive(value_pert, tau_pert, prec);
  _acb_vec_sub(value_pert, value_pert, value, n, prec);
  _acb_vec_scalar_div_arb(value_pert, value_pert, n, eps, prec);
  for (k = 0; k < n; k++)
    {
      acb_set(acb_mat_entry(dth, k, 2), &value_pert[k]);
    }

  _acb_vec_set(th, value, n);
  acb_mat_add_error_mag(dth, error);

  acb_mat_clear(tau_pert);
  arb_clear(eps);
  mag_clear(error);
  _acb_vec_clear(value, n);
  _acb_vec_clear(value_pert, n);
  return res;
}

