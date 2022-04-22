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
  fmpz_poly_struct num_vec[3];
  fmpz_t den;
  slong ell;
  slong nb_roots = 0;
  slong nmax;
  slong* mults;
  fmpq* roots;
  slong k;
  int res;

  I = _fmpz_vec_init(4);
  j = _fmpq_vec_init(3);
  for (k = 0; k < 3; k++) fmpz_poly_init(&num_vec[k]);
  fmpz_init(den);
  /* Init roots later: size depends on ell */
  
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
  roots = _fmpq_vec_init(nmax);
  mults = flint_malloc(nmax * sizeof(slong));
  
  /* Set Igusa invariants */
  igusa_switch_I6prime_fmpz(I, I);
  igusa_from_cov_fmpz(j, I);
  
  /* Evaluate modular equations */
  res = siegel_modeq_eval_Q(num_vec, den, j, ell);
  /* Compute roots */
  if (res) modeq_roots_Q(&nb_roots, roots, mults, &num_vec[0]);
  
  for (k = 0; k < nb_roots; k++)
    {
      res = modeq_isog_invariants_Q(j, num_vec, &roots[k], 3);
      if (!res)
	{
	  flint_printf("FAIL (double root)\n");
	  flint_abort();
	}
      cov_from_igusa_fmpz(I, j);
      igusa_switch_I6_fmpz(I, I);
      
      flint_printf("\n["); fmpz_print(&I[0]); flint_printf(",\n");
      fmpz_print(&I[1]); flint_printf(",\n");
      fmpz_print(&I[2]); flint_printf(",\n");
      fmpz_print(&I[3]); flint_printf("]\n");
    }
  
  _fmpz_vec_clear(I, 4);
  _fmpq_vec_clear(j, 3);
  for (k = 0; k < 3; k++) fmpz_poly_clear(&num_vec[k]);
  fmpz_clear(den);
  _fmpq_vec_clear(roots, nmax);
  flint_free(mults);

  flint_cleanup();  
  return EXIT_SUCCESS;
}
