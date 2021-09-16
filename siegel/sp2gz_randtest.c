
#include "siegel.h"

void
sp2gz_randtest(sp2gz_t m, flint_rand_t state, slong bits)
{
  sp2gz_t n;
  sp2gz_init(n, m->g);

  sp2gz_randtest_triangular(m, state, bits);
  sp2gz_randtest_diagonal(n, state, bits);
  sp2gz_mul(m, m, n);
  sp2gz_J(n);
  sp2gz_mul(m, m, n);
  sp2gz_randtest_triangular(n, state, bits);
  sp2gz_mul(m, m, n);
  sp2gz_J(n);
  sp2gz_mul(m, m, n);
  sp2gz_randtest_diagonal(n, state, bits);
  sp2gz_mul(m, m, n);
  sp2gz_J(n);
  sp2gz_mul(m, m, n);
  sp2gz_randtest_triangular(n, state, bits);
  sp2gz_mul(m, m, n);

  sp2gz_clear(n);
}
