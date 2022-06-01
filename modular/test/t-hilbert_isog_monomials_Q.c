
#include "modular.h"

int main()
{
  slong iter;
  
  flint_printf("hilbert_isog_monomials_Q....");
  fflush(stdout);

  for (iter = 0; iter < 2; iter++)
    {
      fmpz* I1;
      slong ell;
      slong delta;
      slong max_nb_roots = 3;
      slong nb_roots;
      slong wt;
      fmpz* all_M;
      slong nb_M;
      slong* exp_array;
      fmpz_mpoly_ctx_t ctx;
      fmpz_mpoly_t mon;
      fmpz* test_M;
      int print = 1;
      slong k, j, i;
      
      nb_roots = 1;
      wt = 20;
      nb_M = igusa_nb_base_monomials(wt);
      
      I1 = _fmpz_vec_init(4);
      all_M = _fmpz_vec_init(nb_M * max_nb_roots);
      exp_array = flint_malloc(4 * nb_M * sizeof(slong));
      
      fmpz_mpoly_ctx_init(ctx, 4, ORD_LEX);
      fmpz_mpoly_init(mon, ctx);
      test_M = _fmpz_vec_init(nb_M);
      
      if (iter == 0)
	{
	  /* https://beta.lmfdb.org/Genus2Curve/Q/529/a/529/1 */
	  fmpz_set_si(&I1[0], 284);
	  fmpz_set_si(&I1[1], 2401);
	  fmpz_set_si(&I1[2], 246639);
	  fmpz_set_si(&I1[3], -67712);
	  ell = 11;
	  delta = 5;
	}
      else
	{
	  /* https://beta.lmfdb.org/Genus2Curve/Q/841/a/841/1 */
	  fmpz_set_si(&I1[0], 1420);
	  fmpz_set_si(&I1[1], 4201);
	  fmpz_set_si(&I1[2], 1973899);
	  fmpz_set_si(&I1[3], 107648);
	  ell = 7;
	  delta = 8;
	}
      
      igusa_from_IC_fmpz(I1, I1);
      hilbert_isog_monomials_Q(&nb_roots, all_M, &nb_M, exp_array,
			       I1, ell, delta);
      if (nb_roots == 0 || nb_M != igusa_nb_base_monomials(wt))
	{
	  flint_printf("FAIL (roots/weight)\n");
	  flint_printf("nb_roots = %wd, nb_M = %wd\n", nb_roots, nb_M);
	  fflush(stdout);
	  flint_abort();
	}

      if (print)
	{
	  flint_printf("Number of roots: %wd\n", nb_roots);
	  for (k = 0; k < nb_roots; k++)
	    {
	      flint_printf("Value of monomials at root number %wd:\n", k+1);
	      for (j = 0; j < nb_M; j++)
		{		  
		  flint_printf("Exponents");
		  for (i = 0; i < 4; i++) flint_printf(" %wd", exp_array[4*j+i]);
		  flint_printf(": ");
		  fmpz_print(&all_M[nb_M*k + j]);
		  flint_printf("\n");
		}
	    }
	}

      _fmpz_vec_clear(I1, 4);
      _fmpz_vec_clear(all_M, nb_M * max_nb_roots);
      flint_free(exp_array);
      _fmpz_vec_clear(test_M, nb_M);
      fmpz_mpoly_clear(mon, ctx);
      fmpz_mpoly_ctx_clear(ctx);
    }

  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

