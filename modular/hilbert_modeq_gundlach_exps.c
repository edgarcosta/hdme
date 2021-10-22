
#include "modular.h"

void hilbert_modeq_gundlach_exps(slong* e, slong* a, slong* b, slong ell, slong delta)
{
  
  if (delta != 5)
    {
      flint_printf("(hilbert_modeq_gundlach_exps) Error: Gundlach invariants only implemented for discriminant 5\n");
      fflush(stdout);
      flint_abort();
    }
  
  slong wb = 2*10*siegel_nb_cosets(ell);
  slong e, a, b;
  
  e = 2*(wb / 6);

  /* 2e + wb = 10a + 2b, 0\leq b \leq 4 */
  a = (2*e + wb) / 10;
  b = ((2*e + wb) % 10)/2;

  if (2*e + wb != 10*a + 2*b) flint_abort();

  *a_p = a;
  *b_p = b;
  *e_p = e;
}
