
#include "modular.h"

int main()
{
  slong iter;
  
  flint_printf("siegel_2step_isog_monomials_Q....");
  fflush(stdout);

  for (iter = 0; iter < 1; iter++)
    {
      fmpz* I1;
      fmpz* I2;
      slong ell;
      slong max_nb_roots = 2;
      slong nb_roots;
      slong wt;
      fmpz* all_M;
      slong nb_M;
      slong* exp_array;
      fmpz_mpoly_ctx_t ctx;
      fmpz_mpoly_t mon;
      fmpz* test_M;
      slong k;
      int print = 0;
      int res;

      nb_roots = 1;
      wt = 20;
      nb_M = igusa_nb_base_monomials(wt);
      
      I1 = _fmpz_vec_init(4);
      I2 = _fmpz_vec_init(4);
      all_M = _fmpz_vec_init(nb_M * max_nb_roots);
      exp_array = flint_malloc(4 * nb_M * sizeof(slong));
      
      fmpz_mpoly_ctx_init(ctx, 4, ORD_LEX);
      fmpz_mpoly_init(mon, ctx);
      test_M = _fmpz_vec_init(nb_M);
      
      if (iter == 0)
	{
	  /* https://beta.lmfdb.org/Genus2Curve/Q/249/a/249/1,
	     https://beta.lmfdb.org/Genus2Curve/Q/249/a/6723/1 */
	  fmpz_set_si(&I1[0], 108);
	  fmpz_set_si(&I1[1], 57);
	  fmpz_set_si(&I1[2], 2259);
	  fmpz_set_si(&I1[3], -31872);
	  ell = 2;
	  fmpz_set_si(&I2[0], 1932);
	  fmpz_set_si(&I2[1], 87897);
	  fmpz_set_si(&I2[2], 65765571);
	  fmpz_set_si(&I2[3], 860544);
	}
      else
	{
	  /* https://beta.lmfdb.org/Genus2Curve/Q/277/a/277/1,
	     https://beta.lmfdb.org/Genus2Curve/Q/277/a/277/2 */
	  fmpz_set_si(&I1[0], 64);
	  fmpz_set_si(&I1[1], 352);
	  fmpz_set_si(&I1[2], 9552);
	  fmpz_set_si(&I1[3], -1108);
	  ell = 3;
	  fmpz_set_si(&I2[0], 4480);
	  fmpz_set_si(&I2[1], 1370512);
	  fmpz_set_si(&I2[2], 1511819744);
	  fmpz_set_si(&I2[3], -1108);
	}
      
      igusa_from_IC_fmpz(I1, I1);
      igusa_from_IC_fmpz(I2, I2);

      siegel_2step_isog_monomials_Q(&nb_roots, all_M, &nb_M, exp_array,
				    I1, ell);
      if (nb_roots == 0 || nb_M != igusa_nb_base_monomials(wt))
	{
	  flint_printf("FAIL (roots/weight)\n");
	  flint_printf("nb_roots = %wd, nb_M = %wd\n", nb_roots, nb_M);
	  fflush(stdout);
	  flint_abort();
	}

      if (print)
	{
	  flint_printf("Exponents:\n");
	  for (k = 0; k < 4*nb_M; k++)
	    {
	      flint_printf("%wd ", exp_array[k]);
	    }
	  flint_printf("\n");
	}

      for (k = 0; k < nb_M; k++)
	{
	  cov_monomial(mon, &exp_array[4*k], ctx);
	  cov_mpoly_eval_fmpz(&test_M[k], mon, I2, ctx);
	}
      cov_normalize_fmpz_wt1(test_M, test_M, nb_M);

      res = 0;
      for (k = 0; k < nb_roots; k++)
	{
	  if (_fmpz_vec_equal(test_M, &all_M[k*nb_M], nb_M)) res = 1;
	}
      if (!res)
	{
	  flint_printf("FAIL (values)\n");
	  fflush(stdout);
	  flint_abort();
	}      

      _fmpz_vec_clear(I1, 4);
      _fmpz_vec_clear(I2, 4);
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

