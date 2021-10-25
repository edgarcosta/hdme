
#include "modular.h"

void hilbert_modeq_gundlach_exps(slong* e_p, slong* a_p, slong* b_p, slong ell, slong delta)
{
  slong wb = 2*10*hilbert_nb_cosets(ell, delta);
  slong e, a, b;
  
  
  if (delta != 5)
    {
      flint_printf("(hilbert_modeq_gundlach_exps) Error: Gundlach invariants only implemented for discriminant 5\n");
      fflush(stdout);
      flint_abort();
    }
  
  e = 2*(wb / 6);

  /* 2e + wb = 10a + 6b, 0\leq b \leq 4 */
  a = (2*e + wb) / 30;
  b = (2*e + wb) % 30;

  if (b % 6 == 0)
    {
      a = 3*a;
    }
  else if (b % 6 == 2)
    {
      if (b < 20)
	{
	  a = 3*a-1;
	  b = (b+10)/6;
	}
      else
	{
	  a = 3*a+2;
	  b = (b-20)%6;
	}
    }
  else /* b%6 == 4 */
    {
      if (b < 10)
	{
	  a = 3*a-2;
	  b = (b+20)/6;
	}
      else
	{
	  a = 3*a+1;
	  b = (b-10)%6;
	}
    }

  if (2*e + wb != 10*a + 6*b)
    {
      flint_printf("(hilbert_modeq_gundlach_exps) Error: wrong exponents, %wd != 10*%wd + 6*%wd\n", 2*e+wb, a, b);
      fflush(stdout);      
      flint_abort();
    }

  *a_p = a;
  *b_p = b;
  *e_p = e;
}
