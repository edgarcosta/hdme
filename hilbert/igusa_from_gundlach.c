
#include "hilbert.h"

void igusa_from_gundlach(acb_ptr j, acb_srcptr g, slong delta, slong prec)
{
  acb_t g2, f6, f10, den;
  
  char** vars;
  acb_ptr vals;
  fmpq_mpoly_ctx_t ctx;
  fmpq_mpoly_t pol;
  acb_ptr h;
  acb_t j1, j2, j3;

  acb_init(g2);
  acb_init(f6);
  acb_init(f10);
  acb_init(den);
  h = _acb_vec_init(4);
  acb_init(j1);
  acb_init(j2);
  acb_init(j3);
  
  fmpq_mpoly_ctx_init(ctx, 3, ORD_LEX);
  fmpq_mpoly_init(pol, ctx);
  vars = malloc(3 * sizeof(char*));
  vars[0] = malloc(4 * sizeof(char));
  vars[1] = malloc(4 * sizeof(char));
  vars[2] = malloc(4 * sizeof(char));
  vals = _acb_vec_init(3);
  
  if (delta != 5)
    {
      flint_printf("(gundlach_from_igusa) Error: Gundlach invariants only implemented for discriminant 5\n");
      fflush(stdout);
      flint_abort();
    }
  
  /* Set Gundlach covariants up to weighted scalar */
  acb_one(g2);
  acb_inv(f10, &g[0], prec);
  acb_mul(f6, &g[1], f10, prec);

  /* Get polynomials */
  strcpy(vars[0], "g2");
  strcpy(vars[1], "f6");
  strcpy(vars[2], "f10");
  acb_set(&vals[0], g2);
  acb_set(&vals[1], f6);
  acb_set(&vals[2], f10);
  
  gundlach_get_mpoly(pol, (const char**) vars, "I4", delta, ctx);
  fmpq_mpoly_evaluate_all_acb(&h[0], pol, vals, ctx, prec);
  gundlach_get_mpoly(pol, (const char**) vars, "I6", delta, ctx);
  fmpq_mpoly_evaluate_all_acb(&h[1], pol, vals, ctx, prec);
  gundlach_get_mpoly(pol, (const char**) vars, "I10", delta, ctx);
  fmpq_mpoly_evaluate_all_acb(&h[2], pol, vals, ctx, prec);
  gundlach_get_mpoly(pol, (const char**) vars, "I12", delta, ctx);
  fmpq_mpoly_evaluate_all_acb(&h[3], pol, vals, ctx, prec);
  
  /* Set output */
  cov_from_h(h, h, prec);
  igusa_from_cov(j, h, prec);
  
  acb_clear(g2);
  acb_clear(f6);
  acb_clear(f10);
  acb_clear(den);
  _acb_vec_clear(h, 4);
  acb_clear(j1);
  acb_clear(j2);
  acb_clear(j3);
  
  fmpq_mpoly_clear(pol, ctx);
  fmpq_mpoly_ctx_clear(ctx);
  free(vars[0]);
  free(vars[1]);
  free(vars[2]);
  free(vars);  
}
