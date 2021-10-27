
#include "modular.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("modeq_roots_Q....");
  fflush(stdout);

  flint_randinit(state);
  
  for (iter = 0; iter < 1 * arb_test_multiplier(); iter++)
    {
      fmpz_poly_t crv;
      fmpz_poly_t crv_factor;
      fmpz_t coeff;
      slong bits = 5 + n_randint(state, 5);
      fmpq* j;
      fmpz_poly_struct num_vec[3];
      slong nb_roots = 0;
      fmpq* roots;
      fmpz_t den;
      slong ell = 2;
      slong mults[15]; /* 15 = siegel_nb_cosets(2) */
      slong k;

      fmpz_poly_init(crv);
      fmpz_poly_init(crv_factor);
      fmpz_init(coeff);
      j = _fmpq_vec_init(3);
      for (k = 0; k < 3; k++) fmpz_poly_init(&num_vec[k]);
      fmpz_init(den);
      roots = _fmpq_vec_init(15);

      /* Set up curve */
      do
	{
	  fmpz_poly_one(crv);
	  for (k = 0; k < 6; k++)
	    {
	      fmpz_poly_zero(crv_factor);
	      fmpz_randtest(coeff, state, bits);
	      fmpz_poly_set_coeff_fmpz(crv_factor, 0, coeff);
	      fmpz_randtest_not_zero(coeff, state, bits);
	      fmpz_poly_set_coeff_fmpz(crv_factor, 1, coeff);
	      fmpz_poly_mul(crv, crv, crv_factor);
	    }
	}
      while (!fmpz_poly_is_squarefree(crv));
      igusa_from_curve_fmpz(j, crv);

      /* Modular equation of level 2 must have roots (Richelot!) */
      siegel_modeq_eval_Q(num_vec, den, j, ell);
      modeq_roots_Q(&nb_roots, roots, mults, &num_vec[0]);
      if (nb_roots < 15)
	{
	  flint_printf("FAIL (not split)\n");
	  flint_printf("nb_roots = %wd\n", nb_roots);
	  for (k = 0; k < 3; k++)
	    {
	      fmpq_print(&j[k]); flint_printf("\n");
	    }
	  fmpz_poly_print(&num_vec[0]);
	}

      fmpz_poly_clear(crv);
      fmpz_poly_clear(crv_factor);
      fmpz_clear(coeff);
      _fmpq_vec_clear(j, 3);
      for (k = 0; k < 3; k++) fmpz_poly_clear(&num_vec[k]);
      fmpz_clear(den);
      _fmpq_vec_clear(roots, 15);
    }

  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
