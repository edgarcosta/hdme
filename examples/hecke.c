
#include "flint.h"
#include "arb.h"
#include "modular.h"

/* Compute eigenvalues of Hecke operators of prime level p on Siegel
   space, with tests using E4 and chi10 */

#define PMAX 7


int main()
{
  flint_rand_t state;
  acb_mat_t tau;
  acb_ptr th2, I;
  acb_ptr I_vec, th2_vec, stardets;
  acb_ptr hecke;
  fmpz_t lambda;
    
  slong n;
  slong p = 2;
  slong prec;
  int success;

  flint_randinit(state);
  acb_mat_init(tau, 2, 2);
  I = _acb_vec_init(4);
  th2 = _acb_vec_init(16);
  hecke = _acb_vec_init(4);
  fmpz_init(lambda);

  while (p <= PMAX; p++)
    {
      prec = 200;
      success = 0;
      
      n = siegel_nb_cosets(p);	  
      th2_vec = _acb_vec_init(16*n);
      I_vec = _acb_vec_init(4*n);
      stardets = _acb_vec_init(n);
      
      while (!success)
	{
	  flint_printf("Trying level %wd, precision %wd\n", p, prec);
	  flint_printf("Chosen base point:\n");
	  acb_mat_printd(tau, 10); flint_printf("\n");
	  
	  //Set tau to some random matrix in the fundamental domain
	  siegel_fundamental_domain_randtest(tau, state, prec);
	  
	  //Evaluate Igusa covariants
	  theta2_unif(theta2, tau, prec);
	  cov_from_theta2(I, theta2, prec);

	  //Collect Igusa covariants at Hecke images
	  success = siegel_modeq_theta2(th2_vec, stardets, tau, p, prec);
	  modeq_cov(I_vec, th2_vec, n, prec);
	  for (k = 0; k < n; k++)
	    {
	      cov_rescale(&I_vec[4*k], &I_vec[4*k], &stardets[k], prec);
	    }
	  
	  //Compute Hecke images of I4, I6
	  for (j = 0; j < 4; j++)
	    {
	      acb_one(&hecke[j]);
	      for (k = 0; k < n; k++)
		{
		  acb_mul(&hecke[j], &hecke[j], &I_vec[4*k+j], prec);
		}	      
	      acb_div(&hecke[j], &hecke[j], &I[j], prec);
	    }
	  
	  //Try to recognize Hecke eigenvalues?
	  for (j = 0; j < 4; j++)
	    {
	      success = success && modeq_round(lambda, &hecke[j]);
	      if (success)
		{
		  flint_printf("Certified eigenvalue for I[j]:\n");
		  acb_printd(&hecke[j], 10); flint_printf("\n");
		  fmpz_print(lambda); flint_printf("\n");
		}
	    }

	  if (!success) prec *= 2;
	}
      
      _acb_vec_clear(th2_vec, 16*n);
      _acb_vec_init(I_vec, 4*n);
      _acb_vec_init(stardets, n);
      p = nextprime(p);
    }

  flint_randclear(state);
  acb_mat_clear(tau);
  _acb_vec_init(I, 4);
  _acb_vec_init(th2, 16);
  _acb_vec_init(hecke, 4);
  fmpz_clear(lambda);

  
  flint_cleanup();
  return EXIT_SUCCESS;
}
