
#include "hecke.h"

int hecke_all_isog_Q(slong* nb_roots, fmpz* all_I, hecke_t H, fmpz* I, slong prec)
{
  slong k, j, i;
  slong highprec;
  hecke_t H2;
  fmpz* round;
  fmpz_t zero;
  acb_t m, c;
  arf_t rad;
  int res;
  int v = get_hecke_verbose();
  slong weights[4] = IGUSA_HALFWEIGHTS;

  hecke_init(H2, hecke_nb(H));
  round = _fmpz_vec_init(4);
  fmpz_init(zero);
  acb_init(m);
  acb_init(c);
  arf_init(rad);
  
  highprec = hecke_integral_highprec(H, prec);
  if (v) flint_printf("(hecke_all_isog_Q) Chosen high precision: %wd\n", highprec);
  
  res = hecke_set_I_fmpz_with_lowprec(H2, I, hecke_theta2_tau(H), highprec);
  *nb_roots = 0;
  
  for (k = 0; k < hecke_nb(H); k++)
    {
      if (!res) break;
      if (!acb_vec_contains_int(hecke_I_norm(H, k), 4)) continue;

      /* Round to nearest integers; if failure, abort with res = 0 */
      for (j = 0; j < 4; j++)
	{
	  res = res && acb_round(&round[j], rad, &hecke_I_norm(H, k)[j]);
	}
      if (!res) break;

      /* Recompute k-th entry of H at high precision */

      if (v) flint_printf("(hecke_all_isog_Q) Recomputing entry at precision %wd\n", highprec);
      hecke_set_entry(H2, k, hecke_coset(H, k), highprec);
      hecke_normalize_entry(H2, k, I, hecke_norm_ind(H), highprec);

      /* Check norms are zero at low precision */
      for (j = 0; j < 4; j++)
	{
	  acb_sub_fmpz(m, &hecke_I_norm(H2, k)[j], &round[j], highprec);
	  for (i = 0; i < hecke_nb(H); i++)
	    {
	      if (i == k) continue;
	      acb_sub_fmpz(c, &hecke_I_norm(H, i)[j], &round[j], prec);
	      acb_mul(m, m, c, prec);
	    }
	  res = res && acb_round(zero, rad, m);
	  res = res && fmpz_is_zero(zero);
	}      
      if (!res) break;
      
      /* Set result */
      cov_normalize_fmpz(round, round, 4, weights);
      _fmpz_vec_set(&all_I[4*(*nb_roots)], round, 4);
      *nb_roots += 1;
    }
  
  hecke_clear(H2);
  _fmpz_vec_clear(round, 4);
  fmpz_clear(zero);
  acb_clear(m);
  acb_clear(c);
  arf_clear(rad);
  return res;
}
