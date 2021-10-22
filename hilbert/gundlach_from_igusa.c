
#include "hilbert.h"

void gundlach_from_igusa(acb_ptr G, acb_srcptr I, slong delta, slong prec)
{
  acb_t j1, j2, j3, den;
  char** vars;
  acb_ptr vals;
  fmpq_mpoly_ctx_t ctx;
  fmpq_mpoly_t pol;
  acb_t g1, g2;

  acb_init(j1);
  acb_init(j2);
  acb_init(j3);
  acb_init(den);
  acb_init(g1);
  acb_init(g2);

  fmpq_mpoly_ctx_init(ctx, 3, ORD_LEX);
  fmpq_mpoly_init(pol, ctx);
  vars = malloc(3 * sizeof(char*));
  vars[0] = malloc(3 * sizeof(char));
  vars[1] = malloc(3 * sizeof(char));
  vars[2] = malloc(3 * sizeof(char));
  vals = _acb_vec_init(3);
  
  if (delta != 5)
    {
      flint_printf("(gundlach_from_igusa) Error: Gundlach invariants only implemented for discriminant 5\n");
      fflush(stdout);
      flint_abort();
    }

  /* Set variant of Igusa invariants */
  acb_pow_ui(j1, &I[0], 5, prec);
  acb_div(j1, j1, &I[3], prec);
  acb_pow_ui(j2, &I[0], 3, prec);
  acb_mul(j2, j2, &I[1], prec);
  acb_div(j2, j2, &I[3], prec);
  igusa_I6(j3, I, prec);
  acb_mul(j3, j3, &I[0], prec);
  acb_mul(j3, j3, &I[0], prec);
  acb_div(j3, j3, &I[3], prec);

  /* Get polynomials */
  strcpy(vars[0], "a");
  strcpy(vars[1], "b");
  strcpy(vars[2], "c");
  acb_set(&vals[0], j1);
  acb_set(&vals[1], j2);
  acb_set(&vals[2], j3);
  gundlach_get_mpoly(pol, (const char**) vars, "Qnum", delta, ctx);
  fmpq_mpoly_evaluate_all_acb(g2, pol, vals, ctx, prec);
  gundlach_get_mpoly(pol, (const char**) vars, "Qden", delta, ctx);
  fmpq_mpoly_evaluate_all_acb(den, pol, vals, ctx, prec);
  acb_div(g2, g2, den, prec);
  gundlach_get_mpoly(pol, (const char**) vars, "J1num", delta, ctx);
  fmpq_mpoly_evaluate_all_acb(g1, pol, vals, ctx, prec);
  gundlach_get_mpoly(pol, (const char**) vars, "J1den", delta, ctx);
  fmpq_mpoly_evaluate_all_acb(den, pol, vals, ctx, prec);
  acb_div(g1, g1, den, prec);
  acb_mul(g2, g2, g1, prec);  

  /* Set output */
  acb_set(&G[0], g1);
  acb_set(&G[1], g2);
  
  acb_clear(j1);
  acb_clear(j2);
  acb_clear(j3);
  acb_clear(den);
  acb_clear(g1);
  acb_clear(g2);

  fmpq_mpoly_clear(pol, ctx);
  fmpq_mpoly_ctx_clear(ctx);
  free(vars[0]);
  free(vars[1]);
  free(vars[2]);
  free(vars);
  _acb_vec_clear(vals, 3);
}
