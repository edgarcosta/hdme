/*
    Copyright (C) 2021 Jean Kieffer

    This file is part of the hdme library.

    hdme is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License (GPL v3). See
    LICENCE or <http://www.gnu.org/licenses/> for more details.
*/

#include "flint.h"
#include "arb.h"
#include "igusa.h"

/* INPUT: Real and imaginary parts of Igusa-Clebsch invariants I2, I4,
   I6, I10 as well as the desired output precision prec, in the
   following form:

   re2 im2 re4 im4 re6 im6 re10 im10 prec

   where each string can be recognized by either arb_set_str or
   flint_sscanf.

   OUTPUT: Entries \tau_{1,1}, \tau_{2,2}, \tau_{3,3} of a period
   matrix in the Siegel fundamental domain realizing the given
   Igusa-Clebsch invariants up to (weighted) scaling. Output follows
   the format of acb_printd and is separated by newlines.
*/

int main(int argc, char* argv[])
{
  acb_mat_t tau;
  acb_ptr I;
  arb_t re, im;
  acb_t c;
  int res;
  slong k;
  slong prec;
  slong prec10;

  /* Initialize everything */
  acb_mat_init(tau, 2, 2);
  I = _acb_vec_init(4);
  arb_init(re);
  arb_init(im);
  acb_init(c);

  /* Check correct number of arguments */
  if (argc != 10)
    {
      flint_printf("Wrong number of arguments (expected 9)\n");
      flint_abort();
    }

  /* Set up vector of Igusa-Clebsch invariants and prec */
  res = flint_sscanf(argv[9], "%wd", &prec);
  if (!res)
    {
      flint_printf("Could not read integer\n%s\n", argv[9]);
      flint_abort();
    }
  for (k = 0; k < 4; k++)
    {
      res = arb_set_str(re, argv[1+2*k], prec)
	|| arb_set_str(im, argv[1+2*k+1], prec);
      if (res != 0)
	{
	  flint_printf("Could not read real and imaginary parts\n%s\n%s\n",
		       argv[1+2*k], argv[1+2*k]);
	  flint_abort();
	}
      arb_set(acb_realref(c), re);
      arb_set(acb_imagref(c), im);
      acb_set(&I[k], c);
    }

  /* Correct vector of Igusa-Clebsch invariants: program uses I6', not I6 */
  acb_mul(c, &I[0], &I[1], prec);
  acb_addmul_si(c, &I[2], -3, prec);
  acb_div_si(c, c, 2, prec);
  acb_set(&I[2], c);

  /* Compute period matrix */
  tau_from_igusa(tau, I, prec);

  /* Print result in base 10 */
  prec10 = n_clog_2exp(prec, 10);
  acb_printd(acb_mat_entry(tau, 0, 0), prec10); flint_printf("\n");
  acb_printd(acb_mat_entry(tau, 0, 1), prec10); flint_printf("\n");
  acb_printd(acb_mat_entry(tau, 1, 1), prec10); flint_printf("\n");

  /* Clear everything */
  acb_mat_clear(tau);
  _acb_vec_clear(I, 4);
  arb_clear(re);
  arb_clear(im);
  acb_clear(c);  

  flint_cleanup();
  return EXIT_SUCCESS;
}
