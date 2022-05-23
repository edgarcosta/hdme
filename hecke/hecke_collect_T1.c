
#include "hecke.h"

/* See also hecke_collect_siegel */

int hecke_collect_T1(hecke_t H, slong p, slong prec)
{  
  slong k;
  fmpz_mat_t gamma;
  slong nb = hecke_nb(H);
  int res;
  int v = HECKE_VERBOSE;
  
  fmpz_mat_init(gamma, 4, 4);  
  hecke_ell(H) = p;
  hecke_check_nb(H, siegel_nb_T1_cosets(p));
    
  if (v) flint_printf("(hecke_collect_T1) Computing theta constants (%wd)", nb);
  fflush(stdout);

  /* Loop over all cosets to compute desired data */
  for (k = 0; k < nb; k++)
    {
      if (v) hecke_collect_print_status(res, k);
      if (!res) break;
      
      siegel_T1_coset(gamma, k, p);
      res = hecke_set_entry(H, k, gamma, prec);
    }
  if (v) flint_printf("\n");
  
  fmpz_mat_clear(gamma);
  return res;  
}
