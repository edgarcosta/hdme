
#include <assert.h>
#include <stdio.h>
#include "flint.h"
#include "arb.h"
#include "modular.h"

/* Compute eigenvalues of Hecke operators of prime level p on Siegel
   space, with tests using E4 and chi10 */

#define PMAX 29


int main()
{
  flint_rand_t state;
  FILE* file;
  
  acb_mat_t tau;
  acb_ptr th2, I;
  acb_ptr I_vec, th2_vec, stardets;
  acb_ptr hecke;
  fmpz_t term;
  fmpz_t lambda;
  fmpz_mat_t coset;
  fmpz_mat_t block;
  acb_t scal;
    
  slong n;
  slong p = 2;
  slong prec;
  slong j, k;
  int success, r;

  flint_randinit(state);
  acb_mat_init(tau, 2, 2);
  I = _acb_vec_init(4);
  th2 = _acb_vec_init(16);
  hecke = _acb_vec_init(4);
  fmpz_init(term);
  fmpz_init(lambda);
  fmpz_mat_init(coset, 4, 4);
  fmpz_mat_init(block, 2, 2);
  acb_init(scal);

  file = fopen("hecke.txt", "w");

  while (p <= PMAX)
    {
      prec = 100;
      success = 0;
      
      n = siegel_nb_cosets(p);	  
      th2_vec = _acb_vec_init(16*n);
      I_vec = _acb_vec_init(4*n);
      stardets = _acb_vec_init(n);
      
      /* Set tau to some random matrix in the fundamental domain */
      siegel_fundamental_domain_randtest(tau, state, prec);	  
      flint_printf("Chosen base point:\n");
      acb_mat_printd(tau, 10); flint_printf("\n");
            
      while (!success)
	{
	  flint_printf("\nTrying level %wd, precision %wd\n", p, prec);
	  flint_fprintf(file, "%wd ", p);
	  
	  /* Evaluate Igusa covariants */
	  theta2_unif(th2, tau, prec);
	  theta2_renormalize(th2, th2, prec);
	  cov_from_theta2(I, th2, prec);

	  /* Collect Igusa covariants at Hecke images */
	  success = siegel_modeq_theta2(th2_vec, stardets, tau, p, prec);
	  modeq_cov(I_vec, th2_vec, n, prec);
	  for (k = 0; k < n; k++)
	    {
	      acb_inv(&stardets[k], &stardets[k], prec);
	      cov_rescale(&I_vec[4*k], &I_vec[4*k], &stardets[k], prec);
	      
	      /* Adjust depending on Siegel coset */
	      siegel_coset(coset, k, p);
	      fmpz_mat_get_c(block, coset);
	      assert(fmpz_mat_is_zero(block));
	      fmpz_mat_get_d(block, coset);
	      fmpz_mat_det(lambda, block);
	      acb_set_fmpz(scal, lambda);
	      acb_inv(scal, scal, prec);
	      cov_rescale(&I_vec[4*k], &I_vec[4*k], scal, prec);
	    }
	  
	  /* Compute Hecke images of I4, I10 */
	  for (j = 0; j < 4; j++)
	    {
	      acb_zero(&hecke[j]);
	      for (k = 0; k < n; k++)
		{
		  acb_add(&hecke[j], &hecke[j], &I_vec[4*k+j], prec);
		}	      
	      acb_div(&hecke[j], &hecke[j], &I[j], prec);
	    }

	  flint_printf("Hecke\n");
	  for (j = 1; j < 4; j++)
	    {
	      acb_printd(&hecke[j], 10); flint_printf("\n");
	    }
	  acb_set_si(scal, p);
	  acb_pow_si(scal, scal, 2*4-3, prec);
	  acb_mul(&hecke[1], &hecke[1], scal, prec);
	  acb_set_si(scal, p);
	  acb_pow_si(scal, scal, 2*6-3, prec);
	  acb_mul(&hecke[2], &hecke[2], scal, prec);
	  acb_set_si(scal, p);	
	  acb_pow_si(scal, scal, 2*10-3, prec);  
	  acb_mul(&hecke[3], &hecke[3], scal, prec);	  
	  
	  /* Try to recognize Hecke eigenvalues as integers */
	  for (j = 1; j < 4; j++)
	    {
	      r = modeq_round_coeff(lambda, &hecke[j]);
	      success = success && r;
	      if (r)
		{
		  flint_printf("Certified eigenvalue for I[%wd]:\n", j);
		  acb_printd(&hecke[j], 10); flint_printf("\n");
		  fmpz_print(lambda); flint_printf("\n");
		  fmpz_fprint(file, lambda);
		  flint_fprintf(file, " ");
		  if (j == 1) k = 4; else if (j == 2) k = 6; else if (j == 3) k = 10;
		  fmpz_one(lambda);
		  fmpz_set_si(term, p);
		  fmpz_pow_ui(term, term, k-2);
		  fmpz_add(lambda, lambda, term);
		  fmpz_set_si(term, p);
		  fmpz_pow_ui(term, term, k-1);
		  fmpz_add(lambda, lambda, term);
		  fmpz_set_si(term, p);
		  fmpz_pow_ui(term, term, 2*k-3);
		  fmpz_add(lambda, lambda, term);
		  fmpz_print(lambda); flint_printf("\n");
		}
	      
	    }
	  flint_fprintf(file, "\n");

	  if (!success) prec *= 2;
	}
      
      _acb_vec_clear(th2_vec, 16*n);
      _acb_vec_clear(I_vec, 4*n);
      _acb_vec_clear(stardets, n);
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
