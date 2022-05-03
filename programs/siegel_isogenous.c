/*
    Copyright (C) 2022 Jean Kieffer

    This file is part of the hdme library.

    hdme is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License (GPL v3). See
    LICENCE or <http://www.gnu.org/licenses/> for more details.
*/

#include "modular.h"

/* INPUT: Integers I2, I4, I6, I10, ell as strings recognized by
   fmpz_set_str in base 10

   OUTPUT: list of tuples of integers [I2, I4, I6, I10] that are all
   Igusa-Clebsch invariants of rational ell-isogenous abelian
   surfaces, separated by newlines. 
*/

int main(int argc, char* argv[])
{
  fmpz* I;
  fmpq* j;
  slong ell;
  slong nb_roots = 0;
  slong nmax;
  fmpq* all_isog_j;
  slong k;
  int res;

  I = _fmpz_vec_init(4);
  j = _fmpq_vec_init(3);
  /* Init all_isog_j later: size depends on ell */
  
  /* Check correct number of arguments */
  if (argc != 6)
    {
      flint_printf("Wrong number of arguments (expected 5)\n");
      flint_abort();
    }
  
  /* Scan arguments */
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
  nmax = siegel_nb_cosets(ell);
  all_isog_j = _fmpq_vec_init(3*nmax);
  
  /* Set Igusa invariants */
  igusa_switch_I6prime_fmpz(I, I);
  igusa_from_cov_fmpz(j, I);
  
  /* Compute isogenous invariants */
  res = siegel_modeq_isog_invariants_Q(&nb_roots, all_isog_j, j, ell);
  if (!res)
    {
      flint_printf("FAIL\n");
      flint_abort();
    }
  
  for (k = 0; k < nb_roots; k++)
    {
      cov_from_igusa_fmpz(I, &all_isog_j[3*k]);
      igusa_switch_I6_fmpz(I, I);
      
      flint_printf("\n["); fmpz_print(&I[0]); flint_printf(",\n");
      fmpz_print(&I[1]); flint_printf(",\n");
      fmpz_print(&I[2]); flint_printf(",\n");
      fmpz_print(&I[3]); flint_printf("]\n");
    }
  
  _fmpz_vec_clear(I, 4);
  _fmpq_vec_clear(j, 3);
  _fmpq_vec_clear(all_isog_j, 3*nmax);

  flint_cleanup();  
  return EXIT_SUCCESS;
}
