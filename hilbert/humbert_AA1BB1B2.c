
#include "hilbert.h"

void humbert_AA1BB1B2(acb_ptr AA1BB1B2, const acb_t r, const acb_t s, slong delta,
		      slong prec)
{
  fmpq_mpoly_t pol;
  fmpq_mpoly_ctx_t ctx;
  acb_ptr vals;
  char** vars;

  fmpq_mpoly_ctx_init(ctx, 2, ORD_LEX);
  fmpq_mpoly_init(pol, ctx);
  vars = hdme_data_vars_init(2);
  vals = _acb_vec_init(2);
  
  humbert_vars_set(vars, delta);
  acb_set(&vals[0], r);
  acb_set(&vals[1], s);
  
  if (delta == 33)
    {
      humbert_get_mpoly(pol, (const char**) vars, "I2", delta, ctx);
      hdme_data_evaluate_acb(&AA1BB1B2[0], pol, vals, ctx, prec);
      humbert_get_mpoly(pol, (const char**) vars, "I4", delta, ctx);
      hdme_data_evaluate_acb(&AA1BB1B2[1], pol, vals, ctx, prec);
      humbert_get_mpoly(pol, (const char**) vars, "I6", delta, ctx);
      hdme_data_evaluate_acb(&AA1BB1B2[2], pol, vals, ctx, prec);
      humbert_get_mpoly(pol, (const char**) vars, "I10", delta, ctx);
      hdme_data_evaluate_acb(&AA1BB1B2[3], pol, vals, ctx, prec);
      humbert_get_mpoly(pol, (const char**) vars, "den", delta, ctx);
      hdme_data_evaluate_acb(&AA1BB1B2[4], pol, vals, ctx, prec);
    }
  else
    {
      humbert_get_mpoly(pol, (const char**) vars, "A", delta, ctx);
      hdme_data_evaluate_acb(&AA1BB1B2[0], pol, vals, ctx, prec);
      
      /* acb_printd(&vals[0], 30); flint_printf("\n"); */
      /* acb_printd(&vals[1], 30); flint_printf("\n"); */
      /* fmpq_mpoly_print_pretty(pol, (const char**) vars, ctx); flint_printf("\n"); */
      /* acb_printd(&AA1BB1B2[0], 30); flint_printf("\n"); */
      
      humbert_get_mpoly(pol, (const char**) vars, "A1", delta, ctx);
      hdme_data_evaluate_acb(&AA1BB1B2[1], pol, vals, ctx, prec);      
      humbert_get_mpoly(pol, (const char**) vars, "B", delta, ctx);
      hdme_data_evaluate_acb(&AA1BB1B2[2], pol, vals, ctx, prec);      
      humbert_get_mpoly(pol, (const char**) vars, "B1", delta, ctx);
      hdme_data_evaluate_acb(&AA1BB1B2[3], pol, vals, ctx, prec);
      if (delta == 61)
	{   
	  humbert_get_mpoly(pol, (const char**) vars, "B1den", delta, ctx);
	  hdme_data_evaluate_acb(&AA1BB1B2[4], pol, vals, ctx, prec);
	  acb_div(&AA1BB1B2[3], &AA1BB1B2[3], &AA1BB1B2[4], prec);
	}    
      humbert_get_mpoly(pol, (const char**) vars, "B2", delta, ctx);
      hdme_data_evaluate_acb(&AA1BB1B2[4], pol, vals, ctx, prec);
    }
  
  fmpq_mpoly_clear(pol, ctx);
  fmpq_mpoly_ctx_clear(ctx);
  hdme_data_vars_clear(vars, 2);
  _acb_vec_clear(vals, 2);
}
