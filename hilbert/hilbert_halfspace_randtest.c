
#include "hilbert.h"

void hilbert_halfspace_randtest(acb_t t1, acb_t t2, flint_rand_t state, slong prec)
{
  slong mag_bits = 1 + n_randint(state, 5);
  int stop = 0;
  
  arb_randtest_precise(acb_realref(t1), state, prec, mag_bits);
  arb_randtest_precise(acb_realref(t2), state, prec, mag_bits);

  while (!stop)
    {
      arb_randtest_precise(acb_imagref(t1), state, prec, mag_bits);
      arb_randtest_precise(acb_imagref(t2), state, prec, mag_bits);
      stop = arb_is_positive(acb_imagref(t1))
	&& arb_is_positive(acb_imagref(t2));
    }
}
