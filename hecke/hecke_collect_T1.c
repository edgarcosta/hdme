
#include "hecke.h"

/* See also hecke_collect_siegel */

int hecke_collect_T1(hecke_t H, slong ell, slong prec) {
  fmpz_mat_t gamma;
  slong nb = hecke_nb(H);
  int res = 1;
  int v = get_hecke_verbose();

  fmpz_mat_init(gamma, 4, 4);
  hecke_ell(H) = ell;
  hecke_check_nb(H, siegel_nb_T1_cosets(ell));

  if (v) hecke_collect_verbose_start(nb);

  /* Loop over all cosets to compute desired data */
  #pragma omp parallel for shared(res, H)
  for(slong k = 0; k < nb; ++k) {
    if (v) hecke_collect_print_status(res, k, nb);
    if (!res) continue; // OpenMP doesn't allow break, and thus we continue

    siegel_T1_coset(gamma, k, ell);
    int localres = hecke_set_entry(H, k, gamma, prec);
    #pragma omp atomic
    res &= localres;
  }
  if (v) flint_printf("\n");

  hecke_norm_ind(H) = n_pow(ell, 3);
  fmpz_set_si(hecke_norm_all(H), ell);
  fmpz_pow_ui(hecke_norm_all(H), hecke_norm_all(H),
	      3*hecke_nb(H) - (ell*ell + 2*ell + 1));
  hecke_prod_ec(H) = 2*(ell*ell + ell);

  fmpz_mat_clear(gamma);
  return res;
}
