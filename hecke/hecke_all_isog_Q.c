
#include "hecke.h"

int hecke_all_isog_Q(slong* nb_roots, fmpz* all_I, hecke_t H, fmpz* I, slong prec) {
  time_pair start; timestamp_mark(&start);
  slong highprec;
  hecke_t H2;
  int res;
  int v = get_hecke_verbose();
  slong weights[4] = IGUSA_HALFWEIGHTS;
  slong nb = hecke_nb(H);

  hecke_init(H2, nb);

  highprec = hecke_integral_highprec(H, prec);
  if (v) flint_printf("(hecke_all_isog_Q) Chosen high precision: %wd\n", highprec);

  res = hecke_set_I_fmpz_with_lowprec(H2, I, hecke_theta2_tau(H), highprec);
  *nb_roots = 0;

  #pragma omp parallel
	{ // this initializes the private	variables
    fmpz* round;
    fmpz_t zero;
    acb_t m, c;
    arf_t rad;

    round = _fmpz_vec_init(4);
    fmpz_init(zero);
    acb_init(m);
    acb_init(c);
    arf_init(rad);

    // filter out entries that need to be recomputed
    #pragma omp for schedule(dynamic)
    for (slong k = 0; k < nb; ++k) {
      if (!res) continue; // OpenMP doesn't allow break, and thus we continue
      if (!acb_vec_contains_int(hecke_I_norm(H, k), 4)) continue;

      /* Round to nearest integers; if failure, abort with res = 0 */
      int localres = 1;
      for (slong j = 0; j < 4; ++j) {
        localres &= acb_round(&round[j], rad, &hecke_I_norm(H, k)[j]);
      }
#pragma omp atomic
      res &=localres;
      if (!res) continue; // OpenMP doesn't allow break, and thus we continue

      /* Recompute k-th entry of H at high precision */
      if (v) flint_printf("(hecke_all_isog_Q) Recomputing %wd-th entry of H at high precision\n", k);
      time_pair mid; timestamp_mark(&mid);
      hecke_set_entry(H2, k, hecke_coset(H, k), highprec);
      if (v) report("hecke_set_entry", mid);

      timestamp_mark(&mid);
      hecke_normalize_entry(H2, k, I, hecke_norm_ind(H), highprec);
      if (v) report("hecke_normalize_entry", mid);

      /* Check norms are zero at low precision */
      for (slong j = 0; j < 4; ++j) {
        acb_sub_fmpz(m, &hecke_I_norm(H2, k)[j], &round[j], highprec);
        for (slong i = 0; i < nb; ++i) {
          if (i == k) continue;
          acb_sub_fmpz(c, &hecke_I_norm(H, i)[j], &round[j], prec);
          acb_mul(m, m, c, prec);
        }
        localres &= acb_round(zero, rad, m);
        localres &= fmpz_is_zero(zero);
      }
      #pragma omp atomic
      res &=localres;
      if (!res) continue; // OpenMP doesn't allow break, and thus we continue

      cov_normalize_fmpz(round, round, 4, weights);
      #pragma omp critical
      {
        /* Set result */
        _fmpz_vec_set(&all_I[4*(*nb_roots)], round, 4);
        *nb_roots += 1;
      }
    }
    _fmpz_vec_clear(round, 4);
    fmpz_clear(zero);
    acb_clear(m);
    acb_clear(c);
    arf_clear(rad);
  }

  hecke_clear(H2);
  if (v) report_end(start);
  return res;
}
