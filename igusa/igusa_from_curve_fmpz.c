
#include "igusa.h"

void igusa_from_curve_fmpz(fmpz* I, const fmpz_poly_t crv)
{  
  fmpz* ai;
  fmpq_t ai_fmpq[7];
  fmpz_t one;
  fmpq_t ev;
  slong k;
  
  fmpq_mpoly_t pol;
  fmpq_mpoly_ctx_t ctx;
  char** vars;
  fmpz* res;
  
  ai = _fmpz_vec_init(7);
  for (k = 0; k < 7; k++)
    {
      fmpq_init(ai_fmpq[k]);
    }
  
  fmpz_init(one);
  fmpq_init(ev);
  fmpq_mpoly_ctx_init(ctx, 7, ORD_LEX);
  fmpq_mpoly_init(pol, ctx);
  vars = hdme_data_vars_init(7);
  res = _fmpz_vec_init(4);

  hdme_data_vars_set(vars, "a", 0);
  hdme_data_vars_set(vars, "b", 1);
  hdme_data_vars_set(vars, "c", 2);
  hdme_data_vars_set(vars, "d", 3);
  hdme_data_vars_set(vars, "e", 4);
  hdme_data_vars_set(vars, "f", 5);
  hdme_data_vars_set(vars, "g", 6);

  curve_coeffs_fmpz(ai, crv);
  fmpz_one(one);
  for (k = 0; k < 7; k++)
    {
      fmpq_set_fmpz_frac(ai_fmpq[k], &ai[k], one);
    }

  hdme_data_read(pol, (const char**) vars, "igusa/I4", ctx);
  fmpq_mpoly_evaluate_all_fmpq(ev, pol, (fmpq* const*) ai_fmpq, ctx);
  fmpq_numerator(igusa_I4(res), ev);
  hdme_data_read(pol, (const char**) vars, "igusa/I6prime", ctx);
  fmpq_mpoly_evaluate_all_fmpq(ev, pol, (fmpq* const*) ai_fmpq, ctx);
  fmpq_numerator(igusa_I6prime(res), ev);
  hdme_data_read(pol, (const char**) vars, "igusa/I10", ctx);
  fmpq_mpoly_evaluate_all_fmpq(ev, pol, (fmpq* const*) ai_fmpq, ctx);
  fmpq_numerator(igusa_I10(res), ev);
  hdme_data_read(pol, (const char**) vars, "igusa/I2", ctx);
  fmpq_mpoly_evaluate_all_fmpq(ev, pol, (fmpq* const*) ai_fmpq, ctx);
  fmpq_numerator(igusa_I12(res), ev);
  fmpz_mul(igusa_I12(res), igusa_I12(res), igusa_I10(res));

  _fmpz_vec_set(I, res, 4);

  _fmpz_vec_clear(ai, 7);
  for (k = 0; k < 7; k++)
    {
      fmpq_clear(ai_fmpq[k]);
    }
  
  fmpz_clear(one);
  fmpq_clear(ev);
  fmpq_mpoly_clear(pol, ctx);
  fmpq_mpoly_ctx_clear(ctx);
  hdme_data_vars_clear(vars, 7);
  _fmpz_vec_clear(res, 4);  
}
