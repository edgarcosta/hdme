
#include "igusa.h"

void igusa_R2_from_IC_fmpz(fmpq_t R2, fmpz* IC)
{
  fmpq_mpoly_t pol;
  fmpq_mpoly_ctx_t ctx;
  char** vars;
  fmpq* vals;
  slong j;
  
  fmpq_mpoly_ctx_init(ctx, 4, ORD_LEX);
  fmpq_mpoly_init(pol, ctx);
  vars = hdme_data_vars_init(4);
  vals = _fmpq_vec_init(4);

  for (j = 0; j < 4; j++) fmpq_set_fmpz(&vals[j], &IC[j]);

  hdme_data_vars_set(vars, "a", 0);
  hdme_data_vars_set(vars, "b", 1);
  hdme_data_vars_set(vars, "c", 2);
  hdme_data_vars_set(vars, "d", 3);
  
  hdme_data_read(pol, (const char**) vars, "igusa/R2", ctx);
  fmpq_mpoly_evaluate_all_fmpq(res, pol, vals, ctx);

  fmpq_mpoly_clear(pol, ctx);
  fmpq_mpoly_ctx_clear(ctx);
  hdme_data_vars_clear(vars, 4);
  _fmpq_vec_clear(vals, 4);
}
