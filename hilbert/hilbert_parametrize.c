
#include "hilbert.h"

void hilbert_parametrize(acb_ptr I, const acb_t r, const acb_t s, slong delta, slong prec)
{
  acb_t g, h, x, y, z;
  fmpq_mpoly_t pol;
  fmpq_mpoly_ctx_t ctx;
  char** vars;
  acb_ptr vals;

  fmpq_mpoly_ctx_init(ctx, 2, ORD_LEX);
  fmpq_mpoly_init(pol, ctx);
  acb_init(g);
  acb_init(h);
  acb_init(x);
  acb_init(y);
  acb_init(z);
  vars = humbert_vars_init();
  vals = _acb_vec_init(2);

  acb_set(&vals[0], r);
  acb_set(&vals[1], s);

  if (delta == 5)
    {
      strcpy(vars[0], "m");
      strcpy(vars[1], "n");
      hilbert_get_mpoly(pol, (const char**) vars, "g", delta, ctx);
      fmpq_mpoly_evaluate_all_acb(g, pol, vals, ctx, prec);

      strcpy(vars[1], "g");
      acb_set(&vals[1], g);
      
      hilbert_get_mpoly(pol, (const char**) vars, "h", delta, ctx);
      fmpq_mpoly_evaluate_all_acb(h, pol, vals, ctx, prec);      
    }

  else if (delta == 8)
    {      
      strcpy(vars[0], "m");
      strcpy(vars[1], "n");
      hilbert_get_mpoly(pol, (const char**) vars, "rnum", delta, ctx);
      fmpq_mpoly_evaluate_all_acb(g, pol, vals, ctx, prec);      
      hilbert_get_mpoly(pol, (const char**) vars, "rden", delta, ctx);
      fmpq_mpoly_evaluate_all_acb(x, pol, vals, ctx, prec);
      acb_div(g, g, x, prec);

      strcpy(vars[1], "r");
      acb_set(&vals[1], g);
      
      hilbert_get_mpoly(pol, (const char**) vars, "snum", delta, ctx);
      fmpq_mpoly_evaluate_all_acb(h, pol, vals, ctx, prec);      
      hilbert_get_mpoly(pol, (const char**) vars, "sden", delta, ctx);
      fmpq_mpoly_evaluate_all_acb(x, pol, vals, ctx, prec);
      acb_div(h, h, x, prec);
    }

  else if (delta == 12)
    {
      strcpy(vars[0], "f");
      strcpy(vars[1], "g");
      hilbert_get_mpoly(pol, (const char**) vars, "enum", delta, ctx);
      fmpq_mpoly_evaluate_all_acb(g, pol, vals, ctx, prec);      
      hilbert_get_mpoly(pol, (const char**) vars, "eden", delta, ctx);
      fmpq_mpoly_evaluate_all_acb(x, pol, vals, ctx, prec);
      acb_div(g, g, x, prec); /* g is e */

      acb_set(h, &vals[0]);
    }

  else if (delta == 13)
    {    
      strcpy(vars[0], "m");
      strcpy(vars[1], "n");
      hilbert_get_mpoly(pol, (const char**) vars, "g2num", delta, ctx);
      fmpq_mpoly_evaluate_all_acb(x, pol, vals, ctx, prec);      
      hilbert_get_mpoly(pol, (const char**) vars, "g2den", delta, ctx);
      fmpq_mpoly_evaluate_all_acb(y, pol, vals, ctx, prec);
      acb_div(x, x, y, prec); /* x is g2 */
      acb_div_si(y, x, 54, prec); /* y is g1 */
      acb_mul(z, &vals[0], y, prec); /* z is h1 */
      acb_set_si(h, 64);
      acb_div_si(h, h, 27, prec);
      acb_add(h, h, z, prec);
      acb_set_si(g, -1);
      acb_div_si(g, g, 54, prec);
      acb_add(g, g, y, prec);
    }

  else if (delta == 17)
    {
      strcpy(vars[0], "m");
      strcpy(vars[1], "n");
      hilbert_get_mpoly(pol, (const char**) vars, "gnum", delta, ctx);
      fmpq_mpoly_evaluate_all_acb(g, pol, vals, ctx, prec);      
      hilbert_get_mpoly(pol, (const char**) vars, "gden", delta, ctx);
      fmpq_mpoly_evaluate_all_acb(x, pol, vals, ctx, prec);
      acb_div(g, g, x, prec);

      strcpy(vars[1], "g");
      acb_set(&vals[1], g);

      hilbert_get_mpoly(pol, (const char**) vars, "h", delta, ctx);
      fmpq_mpoly_evaluate_all_acb(h, pol, vals, ctx, prec); 
    }

  else
    {
      flint_printf("(hilbert_parametrize) Error: no Hilbert parametrization for discriminant %wd\n", delta);
      fflush(stdout);
      flint_abort();
    }

  humbert_parametrize(I, g, h, delta, prec);

  fmpq_mpoly_clear(pol, ctx);
  fmpq_mpoly_ctx_clear(ctx);
  acb_clear(g);
  acb_clear(h);
  acb_clear(x);
  acb_clear(y);
  acb_clear(z);
  humbert_vars_clear(vars);
  _acb_vec_clear(vals, 2);
}
