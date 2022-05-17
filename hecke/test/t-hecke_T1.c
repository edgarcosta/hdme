
#include "hecke.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("hecke_T1....");
  fflush(stdout);
  flint_randinit(state);

  for (iter = 0; iter < 10 * arb_test_multiplier(); iter++)
    {      
      slong ell = 2;
      slong ellmax = 4;
      acb_mat_t tau;
      hecke_t H;
      acb_ptr val;
      acb_t r;
      fmpz* eig;
      slong k;
      slong nb;
      slong prec = 500 + n_randint(state, 1000);
      
      acb_mat_init(tau, 2, 2);
      acb_init(r);
      eig = _fmpz_vec_init(4);
      
      siegel_fundamental_domain_randtest(tau, state, prec);
      
      while (ell < ellmax)
	{
	  nb = siegel_nb_T1_cosets(ell);
	  hecke_init(H, nb);
	  val = _acb_vec_init(nb);

	  hecke_set_T1(H, tau, ell, prec);
	  
	  /* Test expected eigenvalues for I4, I6' */
	  for (k = 0; k < nb; k++) acb_set(&val[k], &hecke_I(H, k)[0]); /* I4 */
	  hecke_operator(r, H, val, ell, 4, 0, prec);
	  acb_div(r, r, &hecke_I_tau(H)[0], prec);
	  
	  hecke_eigenvalues_eisenstein_p2(eig, 4, ell);
	  if (!acb_contains_fmpz(r, &eig[1]))
	    {
	      flint_printf("FAIL (I4 eigenvalue)\n");
	      flint_printf("ell = %wd\n", ell);
	      acb_printd(r, 30); flint_printf("\n");
	      fmpz_print(&eig[1]); flint_printf("\n");
	      fflush(stdout);
	      flint_abort();
	    }
	  
	  for (k = 0; k < nb; k++) acb_set(&val[k], &hecke_I(H, k)[1]); /* I6' */
	  hecke_operator(r, H, val, ell, 6, 0, prec);
	  acb_div(r, r, &hecke_I_tau(H)[1], prec);
	  
	  hecke_eigenvalues_eisenstein_p2(eig, 6, ell);
	  if (!acb_contains_fmpz(r, &eig[1]))
	    {
	      flint_printf("FAIL (I6' eigenvalue)\n");
	      flint_printf("ell = %wd\n", ell);
	      acb_printd(r, 30); flint_printf("\n");
	      fmpz_print(&eig[1]); flint_printf("\n");
	      fflush(stdout);
	      flint_abort();
	    }
	  
	  hecke_clear(H);
	  _acb_vec_clear(val, nb);
	  ell = n_nextprime(ell, 1);
	}

      acb_mat_clear(tau);
      acb_clear(r);
      _fmpz_vec_clear(eig, 4);
    }
  
  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

