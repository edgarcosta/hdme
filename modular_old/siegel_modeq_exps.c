
#include "modular.h"

void siegel_modeq_exps(slong* e_p, slong* a_p, slong* b_p, slong ell)
{
  slong wl = 20 * siegel_nb_cosets(ell);
  slong e, a, b;
  
  e = wl / 6;

  /* 4e + wl = 10a + 4b, 0\leq b \leq 4 */
  a = (4*e + wl) / 20;
  b = (4*e + wl) % 20;
  if (b%4 == 0)
    {
      a = 2*a;
      b = b/4;
    }
  else if (b < 10) /* b is 2 or 6 */
    {
      b = (b + 10)/4;
      a = 2*a - 1;
    }
  else /* b is 10 or 14 or 18 */
    {
      b = (b - 10)/4;
      a = 2*a + 1;
    }
  if (4*e + wl != 10*a + 4*b) flint_abort();

  *a_p = a;
  *b_p = b;
  *e_p = e;
}
