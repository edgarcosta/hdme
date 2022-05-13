
#include "hecke.h"

int hecke_set_siegel(hecke_t H, const acb_mat_t tau, slong ell, slong prec)
{  
  slong k;
  fmpz_mat_t gamma;
  slong nb = hecke_nb(H);
  int res;
  int v = HECKE_VERBOSE;
  
  fmpz_mat_init(gamma, 4, 4);
  
  hecke_ell(H) = ell;
  acb_mat_set(hecke_tau(H), tau);
  
  /* Sanity check: number of initialized elements */
  if (nb != siegel_nb_cosets(ell))
    {
      flint_printf("(hecke_set_siegel) Error: Hecke data structure initialized with %wd slots instead of the expected %wd\n", nb, siegel_nb_cosets(ell));
      fflush(stdout);
      flint_abort();
    }
  
  if (v) flint_printf("(hecke_set_siegel) Computing theta constants (%wd)", nb);
  fflush(stdout);

  /* Loop over all cosets to compute desired data */
  for (k = 0; k < nb; k++)
    {
      /* Print some status info? */
      if (v)
	{
	  if ((k+1) % 100 == 0)
	    {
	      flint_printf("\n(hecke_set_siegel) (%wd/%wd)", k+1, nb);
	    }
	  flint_printf("."); fflush(stdout);
	}

      /* Compute isogenous period matrices */
      siegel_coset(gamma, k, ell);
      res = hecke_set_entry(H, k, gamma, prec);
      
      if (!res)
	{
	  flint_printf("(hecke_set_siegel) Warning: computation aborted due to low precision\n");
	  break;
	}
    }     
  if (v) flint_printf("\n");

  fmpz_mat_clear(gamma);
  return res;  
}
