/*
    Copyright (C) 2021 Jean Kieffer

    This file is part of the hdme library.

    hdme is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License (GPL v3). See
    LICENCE or <http://www.gnu.org/licenses/> for more details.
*/

#include "flint.h"
#include "arb.h"
#include "modular.h"

/* INPUT: Integers I2, I4, I6, I10, ell as strings recognized by
   fmpz_set_str in base 10.

   OUTPUT: Three polynomials in ZZ[X] (numerators) and one integer
   (denominator) giving the modular equations of Siegel type of level
   ell evaluated at that point. Output format follows
   fmpz_poly_print_pretty and fmpz_print, and is separated by
   newlines. 

   Recall that the Igusa invariants used are
   j1 = I4*I6'/I10, j2 = I2*I4^2/I10, j3 = I4^5/I10^2.
*/

int main(int argc, char* argv[])
{
  fmpz_poly_struct num_vec[3];
  fmpz_t den;
  fmpz* I;
  fmpq* j;
  slong ell = 0;
  slong k;
  int res;

  /* Initialize everything */
  for (k = 0; k < 3; k++) fmpz_poly_init(&num_vec[k]);
  fmpz_init(den);
  I = _fmpz_vec_init(4);
  j = _fmpq_vec_init(3);
  
  /* Check correct number of arguments */
  if (argc != 6)
    {
      flint_printf("Wrong number of arguments (expected 5)\n");
      flint_abort();
    }
  
  /* Set entries of tau and prec */
  res = flint_sscanf(argv[5], "%wd", &ell);
  if (!res || !n_is_prime(ell) || ell > 100)
    {
      flint_printf("Not a small prime:\n%s\n", argv[5]);
      flint_abort();
    }  
  for (k = 0; k < 4; k++)
    {
      res = fmpz_set_str(&I[k], argv[1+k], 10);
      if (res != 0)
	{
	  flint_printf("Could not read integer\n%s\n",
		       argv[1+k]);
	  flint_abort();
	}
    }

  /* Set Igusa invariants */
  igusa_from_cov_fmpz(j, I); /* This will fail if I10 is zero */
  /* Evaluate modular equations */
  siegel_modeq_eval_Q(num_vec, den, j, ell);
  /* Print result (huge!) */
  for (k = 0; k < 3; k++)
    {
      fmpz_poly_print_pretty(&num_vec[k], "x"); flint_printf("\n");
    }
  fmpz_print(den); flint_printf("\n");
  
  /* Clear everything */
  for (k = 0; k < 3; k++) fmpz_poly_clear(&num_vec[k]);
  fmpz_clear(den);
  _fmpz_vec_clear(I, 4);
  _fmpq_vec_clear(j, 3);
  
  flint_cleanup();
  return EXIT_SUCCESS;
}

