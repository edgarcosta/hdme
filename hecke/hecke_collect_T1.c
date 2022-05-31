
#include "hecke.h"

/* See also hecke_collect_siegel */

int hecke_collect_T1(hecke_t H, slong p, slong prec)
{  
  slong k;
  fmpz_mat_t gamma;
  slong nb = hecke_nb(H);
  int res = 1;
  int v = HECKE_VERBOSE;
  
  fmpz_mat_init(gamma, 4, 4);  
  hecke_ell(H) = p;
  hecke_check_nb(H, siegel_nb_T1_cosets(p));
    
  if (v) hecke_collect_verbose_start(nb);

  /* Loop over all cosets to compute desired data */
  for (k = 0; k < nb; k++)
    {
      if (v) hecke_collect_print_status(res, k, nb);
      if (!res) break;
      
      siegel_T1_coset(gamma, k, p);
      res = hecke_set_entry(H, k, gamma, prec);
    }
  if (v) flint_printf("\n");

  fmpz_set_si(hecke_normalize(H), p);
  fmpz_pow_ui(hecke_normalize(H), hecke_normalize(H), 3*hecke_nb(H) - (p*p + 2*p + 1));
  hecke_prod_ec(H) = 2*(p*p + p);
  
  fmpz_mat_clear(gamma);
  return res;  
}
