
#include "igusa.h"

void igusa_scalar_covariants_fmpz(fmpz* I, const fmpz_poly_t crv)
{  
  fmpz* ai;
  fmpq** ai_fmpq;
  fmpz_t one;
  fmpq_t ev;
  slong k;
  
  fmpq_mpoly_t pol;
  fmpq_mpoly_ctx_t ctx;
  char** vars;
  fmpz* res;
  
  ai = _fmpz_vec_init(7);
  ai_fmpq = flint_malloc(7 * sizeof(fmpq*));
  for (k = 0; k < 7; k++)
    {
      ai_fmpq[k] = flint_malloc(sizeof(fmpq));
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

  hdme_data_read(pol, (const char**) vars, "igusa/I2", ctx);
  
  /* fmpz_poly_print(crv);
  fmpq_print(ev); flint_printf("\n");
  fmpq_mpoly_print_pretty(pol, (const char**) vars, ctx); flint_printf("\n");
  for (k = 0; k < 7; k++)
  {fmpq_print(ai_fmpq[k]); flint_printf("\n");} */
  
  fmpq_mpoly_evaluate_all_fmpq(ev, pol, (fmpq* const*) ai_fmpq, ctx);
  fmpq_numerator(&res[0], ev);
  hdme_data_read(pol, (const char**) vars, "igusa/I4", ctx);
  fmpq_mpoly_evaluate_all_fmpq(ev, pol, (fmpq* const*) ai_fmpq, ctx);
  fmpq_numerator(&res[1], ev);
  hdme_data_read(pol, (const char**) vars, "igusa/I6prime", ctx);
  fmpq_mpoly_evaluate_all_fmpq(ev, pol, (fmpq* const*) ai_fmpq, ctx);
  fmpq_numerator(&res[2], ev);
  hdme_data_read(pol, (const char**) vars, "igusa/I10", ctx);
  fmpq_mpoly_evaluate_all_fmpq(ev, pol, (fmpq* const*) ai_fmpq, ctx);
  fmpq_numerator(&res[3], ev);

  _fmpz_vec_set(I, res, 4);

  _fmpz_vec_clear(ai, 7);
  for (k = 0; k < 7; k++)
    {
      fmpq_clear(ai_fmpq[k]);
      flint_free(ai_fmpq[k]);
    }
  flint_free(ai_fmpq);
  
  fmpz_clear(one);
  fmpq_clear(ev);
  fmpq_mpoly_clear(pol, ctx);
  fmpq_mpoly_ctx_clear(ctx);
  hdme_data_vars_clear(vars, 7);
  _fmpz_vec_clear(res, 4);  
}
