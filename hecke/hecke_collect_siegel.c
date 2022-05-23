
#include "hecke.h"

int hecke_collect_siegel(hecke_t H, slong ell, slong prec)
{  
  slong k;
  fmpz_mat_t gamma;
  slong nb = hecke_nb(H);
  int res = 1;
  int v = HECKE_VERBOSE;
  
  fmpz_mat_init(gamma, 4, 4);
  
  hecke_ell(H) = ell;
  hecke_check_nb(H, siegel_nb_cosets(ell));
  
  if (v) flint_printf("(hecke_collect_siegel) Computing theta constants (%wd)", nb);
  fflush(stdout);

  /* Loop over all cosets to compute desired data */
  for (k = 0; k < nb; k++)
    {
      if (v) hecke_collect_print_status(res, k);
      if (!res) break;
      
      siegel_coset(gamma, k, ell);
      res = hecke_set_entry(H, k, gamma, prec);
    }     
  if (v) flint_printf("\n");

  fmpz_mat_clear(gamma);
  return res;  
}
