
#include "hecke.h"

int hecke_set_T1(hecke_t H, const acb_mat_t tau, slong p, slong prec)
{  
  slong k;
  fmpz_mat_t gamma;
  slong nb = hecke_nb(H);
  int res;
  int v = HECKE_VERBOSE;
  
  fmpz_mat_init(gamma, 4, 4);  
  hecke_ell(H) = p;
  
  /* Sanity check: number of initialized elements */
  if (nb != siegel_nb_T1_cosets(p))
    {
      flint_printf("(hecke_set_T1) Error: Hecke data structure initialized with %wd slots instead of the expected %wd\n", nb, siegel_nb_T1_cosets(p));
      fflush(stdout);
      flint_abort();
    }

  res = hecke_set_tau(H, tau, prec);
  
  if (v) flint_printf("(hecke_set_T1) Computing theta constants (%wd)", nb);
  fflush(stdout);

  /* Loop over all cosets to compute desired data */
  for (k = 0; k < nb; k++)
    {
      /* Check res */
      if (!res)
	{
	  flint_printf("(hecke_set_T1) Warning: computation aborted due to low precision\n");
	  break;
	}
      /* Print some status info? */
      if (v)
	{
	  if ((k+1) % 100 == 0)
	    {
	      flint_printf("\n(hecke_set_T1) (%wd/%wd)", k+1, nb);
	    }
	  flint_printf("."); fflush(stdout);
	}
      /* Compute isogenous period matrices */
      siegel_T1_coset(gamma, k, p);
      res = hecke_set_entry(H, k, gamma, prec);
    }     
  if (v) flint_printf("\n");
  
  fmpz_mat_clear(gamma);
  return res;  
}
