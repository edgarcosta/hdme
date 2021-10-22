
#include "hilbert.h"

void gundlach_from_cov(acb_ptr g, acb_srcptr G, slong delta, slong prec)
{
  acb_t g1, g2;
  
  if (delta != 5)
    {
      flint_printf("(gundlach_from_cov) Error: Gundlach invariants only implemented for discriminant 5\n");
      fflush(stdout);
      flint_abort();
    }
  
  acb_init(g1);
  acb_init(g2);

  acb_pow_ui(g1, &G[0], 5, prec);
  acb_div(g1, g1, &G[2], prec);
  acb_pow_ui(g2, &G[0], 2, prec);
  acb_mul(g2, g2, &G[1], prec);
  acb_div(g2, g2, &G[2], prec);

  acb_set(&g[0], g1);
  acb_set(&g[1], g2);

  acb_clear(g1);
  acb_clear(g2);
}
