
#include "igusa.h"

void igusa_from_curve(acb_ptr I, const acb_poly_t crv, slong prec)
{
  acb_ptr ai;
  fmpq_mpoly_t pol;
  fmpq_mpoly_ctx_t ctx;
  char** vars;
  acb_ptr res;
  
  ai = _acb_vec_init(7);
  fmpq_mpoly_ctx_init(ctx, 7, ORD_LEX);
  fmpq_mpoly_init(pol, ctx);
  vars = hdme_data_vars_init(7);
  res = _acb_vec_init(4);

  hdme_data_vars_set(vars, "a", 0);
  hdme_data_vars_set(vars, "b", 1);
  hdme_data_vars_set(vars, "c", 2);
  hdme_data_vars_set(vars, "d", 3);
  hdme_data_vars_set(vars, "e", 4);
  hdme_data_vars_set(vars, "f", 5);
  hdme_data_vars_set(vars, "g", 6);

  curve_coeffs(ai, crv);

  hdme_data_read(pol, (const char**) vars, "igusa/I4", ctx);
  hdme_data_evaluate_acb(igusa_I4(res), pol, ai, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "igusa/I6prime", ctx);
  hdme_data_evaluate_acb(igusa_I6prime(res), pol, ai, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "igusa/I10", ctx);
  hdme_data_evaluate_acb(igusa_I10(res), pol, ai, ctx, prec);
  hdme_data_read(pol, (const char**) vars, "igusa/I2", ctx);
  hdme_data_evaluate_acb(igusa_I12(res), pol, ai, ctx, prec);
  acb_mul(igusa_I12(res), igusa_I12(res), igusa_I10(res), prec);

  _acb_vec_set(I, res, 4);

  _acb_vec_clear(ai, 7);
  fmpq_mpoly_clear(pol, ctx);
  fmpq_mpoly_ctx_clear(ctx);
  hdme_data_vars_clear(vars, 7);
  _acb_vec_clear(res, 4);  
}
