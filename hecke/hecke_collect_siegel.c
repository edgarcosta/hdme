
#include "hecke.h"

int hecke_collect_siegel(hecke_t H, slong ell, slong prec)
{  
  slong k;
  fmpz_mat_t gamma;
  slong nb = hecke_nb(H);
  int res = 1;
  int v = get_hecke_verbose();
  
  fmpz_mat_init(gamma, 4, 4);
  
  hecke_ell(H) = ell;
  hecke_check_nb(H, siegel_nb_cosets(ell));
  
  if (v) hecke_collect_verbose_start(nb);

  /* Loop over all cosets to compute desired data */
  for (k = 0; k < nb; k++)
    {
      if (v) hecke_collect_print_status(res, k, nb);
      if (!res) break;
      
      siegel_coset(gamma, k, ell);
      res = hecke_set_entry(H, k, gamma, prec);
    }     
  if (v) flint_printf("\n");

  /* Set normalization factor */
  hecke_norm_ind(H) = n_pow(ell, 2);
  fmpz_set_si(hecke_norm_all(H), ell);
  fmpz_pow_ui(hecke_norm_all(H), hecke_norm_all(H), 2*hecke_nb(H) - (ell*ell + ell + 2));
  hecke_prod_ec(H) = n_pow(ell+1, 2);
  
  fmpz_mat_clear(gamma);
  return res;  
}
