
#include "modular.h"

/* Test: compute the (9,3,3)-isogeny Shiva discovered in the LMFDB */
int main()
{  
  fmpz* I;
  fmpq* j;
  slong ell = 3;
  slong nb_roots = 0;
  fmpq* all_isog_j;
  slong k;
  
  flint_printf("siegel_modeq_2step....");
  fflush(stdout);

  /* Should be only one isogenous curve */
  I = _fmpz_vec_init(4);
  j = _fmpq_vec_init(3);
  all_isog_j = _fmpq_vec_init(3);
  
  fmpz_set_si(&I[0], 64);
  fmpz_set_si(&I[1], 352);
  fmpz_set_si(&I[2], 9552);
  fmpz_set_si(&I[3], -1108);

  /* We need I6', not I6 */
  igusa_I6prime_fmpz(&I[2], I);
  igusa_from_cov_fmpz(j, I);

  siegel_modeq_2step_isog_invariants_Q(&nb_roots, all_isog_j, j, ell);
  
  /* Web page is https://beta.lmfdb.org/Genus2Curve/Q/277/a/277/2 */
  fmpz_set_si(&I[0], 4480);
  fmpz_set_si(&I[1], 1370512);
  fmpz_set_si(&I[2], 1511819744);
  fmpz_set_si(&I[3], -1108);
  /* We need I6', not I6 */
  igusa_I6prime_fmpz(&I[2], I);
  igusa_from_cov_fmpz(j, I);
  
  if (nb_roots != 1 ||
      !fmpq_equal(&j[0], &all_isog_j[0]) ||
      !fmpq_equal(&j[1], &all_isog_j[1]) ||
      !fmpq_equal(&j[2], &all_isog_j[2]))
    {      
      flint_printf("FAIL\n");
      flint_printf("nb_roots: %wd\n", nb_roots);
      for (k = 0; k < 3*nb_roots; k++)
	{
	  fmpq_print(&all_isog_j[k]); flint_printf("\n");
	}
      fflush(stdout);
      flint_abort();
    }
  
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
      

      
      
