
#include <assert.h>
#include <stdio.h>
#include <flint/flint.h>
#include "arb.h"
#include "hecke.h"

/* Compute eigenvalues of Hecke operators of prime level p on Siegel
   space, with tests using E4 and chi10 */

#define PMAX 29

int main()
{
  flint_rand_t state;
  FILE* file;
  
  acb_mat_t tau;
  hecke_t H;
  acb_ptr val;
  acb_ptr hecke;
  fmpz_t lambda;
  fmpz_t test;
    
  slong n;
  slong p = 2;
  slong prec;
  slong j, k;
  int success, r;
  slong weights[4] = {4, 6, 10, 12};
  slong wt;

  flint_randinit(state);
  acb_mat_init(tau, 2, 2);
  hecke = _acb_vec_init(4);
  fmpz_init(lambda);
  fmpz_init(test);

  file = fopen("hecke.txt", "w");

  while (p <= PMAX)
    {
      prec = 100;
      success = 0;      
      n = siegel_nb_cosets(p);

      hecke_init(H, n);
      val = _acb_vec_init(n);
      
      /* Set tau to some random matrix in the fundamental domain */
      siegel_fundamental_domain_randtest(tau, state, prec);	  
      flint_printf("Chosen base point:\n");
      acb_mat_printd(tau, 10); flint_printf("\n");
            
      while (!success)
	{
	  flint_printf("\nTrying level %wd, precision %wd\n", p, prec);
	  flint_fprintf(file, "%wd ", p);
	  
	  /* Evaluate Igusa covariants */
	  success = hecke_set_siegel(H, tau, p, prec);
	  
	  /* Compute Hecke images of I4, I6' I10, I12 */
	  for (j = 0; j < 4; j++)
	    {
	      wt = weights[j];
	      for (k = 0; k < n; k++) acb_set(&val[k], &hecke_I(H, k)[j]);
	      hecke_operator(&hecke[j], H, val, p, wt, 0, prec);
	      acb_div(&hecke[j], &hecke[j], &hecke_I_tau(H)[j], prec);
	      
	      flint_printf("Hecke (weight %wd)\n", wt);
	      acb_printd(&hecke[j], 10); flint_printf("\n");
	      
	      r = modeq_round_coeff(lambda, &hecke[j]);
	      success = success && r;
	      if (r)
		{
		  flint_printf("Certified eigenvalue for I%wd:\n", wt);
		  fmpz_print(lambda); flint_printf("\n");
		  fmpz_fprint(file, lambda);
		  hecke_eigenvalue_eisenstein_p(test, wt, p);
		  assert (fmpz_equal(lambda, test));
		}	      
	    }
	  flint_fprintf(file, "\n");
	  
	  if (!success) prec *= 2;
	}

      hecke_clear(H);
      _acb_vec_clear(val, n);
      p = n_nextprime(p, 1);
    }

  fclose(file);

  flint_randclear(state);
  acb_mat_clear(tau);
  _acb_vec_clear(I, 4);
  _acb_vec_clear(th2, 16);
  _acb_vec_clear(hecke, 4);
  fmpz_clear(term);
  fmpz_clear(lambda);
  fmpz_mat_clear(coset);
  fmpz_mat_clear(block);
  acb_clear(scal);
  
  flint_cleanup();
  return EXIT_SUCCESS;
}
