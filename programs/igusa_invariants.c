/*
    Copyright (C) 2022 Jean Kieffer

    This file is part of the hdme library.

    hdme is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License (GPL v3). See
    LICENCE or <http://www.gnu.org/licenses/> for more details.
*/

#include "../igusa.h"

/* INPUT: Entries \tau_{1,1}, \tau_{2,2}, \tau_{3,3} of a small period
   matrix, as well as the desired output precision prec, in the
   following form:
   
   re11 im11 re12 im12 re13 im13 prec

   where each string can be recognized by either arb_set_str or
   flint_sscanf. The matrix is assumed to lie in the Siegel
   fundamental domain. (If necessary, use the function
   siegel_fundamental_domain in siegel.h to reduce the input).

   OUTPUT: Values of Igusa-Clebsch covariants

   I4, I6' = (I2*I4 - 3*I6)/2, I10, I12=I2*I10
   
   at \tau. (These combinations are modular forms, and hence are
   defined everywhere.) Output follows the format of acb_printn and is
   separated by newlines.
*/

int main(int argc, char* argv[])
{
  acb_mat_t tau;
  acb_ptr tau_entries;
  acb_ptr theta2;
  acb_ptr I;
  arb_t re, im;
  acb_t c;
  arb_t tol;
  slong prec, prec10;
  slong k;
  int res;

  /* Initialize everything */
  acb_mat_init(tau, 2, 2);
  I = _acb_vec_init(4);
  tau_entries = _acb_vec_init(3);
  theta2 = _acb_vec_init(16);
  arb_init(re);
  arb_init(im);
  acb_init(c);
  arb_init(tol);

  /* Check correct number of arguments */
  if (argc != 8)
    {
      flint_printf("Wrong number of arguments (expected 7)\n");
      flint_abort();
    }
  
  /* Set entries of tau and prec */
  res = flint_sscanf(argv[7], "%wd", &prec);
  if (!res)
    {
      flint_printf("Could not read integer\n%s\n", argv[7]);
      flint_abort();
    }
  for (k = 0; k < 3; k++)
    {
      res = arb_set_str(re, argv[1+2*k], prec)
	|| arb_set_str(im, argv[2+2*k], prec);
      if (res != 0)
	{
	  flint_printf("Could not read real and imaginary parts\n%s\n%s\n",
		       argv[1+2*k], argv[1+2*k]);
	  flint_abort();
	}
      arb_set(acb_realref(c), re);
      arb_set(acb_imagref(c), im);
      acb_set(&tau_entries[k], c);
    }
  acb_set(acb_mat_entry(tau, 0, 0), &tau_entries[0]);
  acb_set(acb_mat_entry(tau, 0, 1), &tau_entries[1]);
  acb_set(acb_mat_entry(tau, 1, 0), &tau_entries[1]);
  acb_set(acb_mat_entry(tau, 1, 1), &tau_entries[2]);

  /* Check tau is in the fundamental domain, with some tolerance at
     the boundary */
  arb_one(tol);
  arb_mul_2exp_si(tol, tol, -prec/4);
  res = siegel_is_in_fundamental_domain(tau, tol, prec);
  if (!res)
    {
      flint_printf("Not in Siegel fundamental domain:\n");
      acb_mat_printd(tau, 10);
      flint_printf("\n");
      flint_abort();
    }

  /* Compute projective vector of theta constants */
  theta2_unif(theta2, tau, prec);
  /* Renormalize by correct scalar factor */
  theta2_renormalize(theta2, theta2, prec);
  /* Compute associated Igusa-Clebsch invariants */
  igusa_from_theta2(I, theta2, prec);
  igusa_streng(I, I, prec);
  
  /* Print result in base 10 */
  prec10 = n_clog_2exp(prec, 10);
  for (k = 0; k < 4; k++)
    {
      acb_printd(&I[k], prec10); flint_printf("\n");
    }
  
  /* Clear everything */
  acb_mat_clear(tau);
  _acb_vec_clear(I, 4);
  _acb_vec_clear(tau_entries, 3);
  _acb_vec_clear(theta2, 16);
  arb_clear(re);
  arb_clear(im);
  acb_clear(c);
  arb_clear(tol);
  
  flint_cleanup();
  return EXIT_SUCCESS;
}
