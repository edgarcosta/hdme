
#include "hecke.h"

int hecke_set_hilbert(hecke_t H, acb_srcptr t, const fmpz_poly_t beta,
		      slong ell, slong delta, slong prec)
{  
  slong k;
  fmpz_poly_mat_t m;
  fmpz_mat_t gamma;
  slong nb = hecke_nb(H);
  int res;
  int v = HECKE_VERBOSE;

  fmpz_poly_mat_init(m, 2, 2);
  fmpz_mat_init(gamma, 4, 4);

  /* Set Hecke context */
  hecke_ell(H) = ell;
  fmpz_poly_set(hecke_beta(H), beta);
  _acb_vec_set(hecke_t1t2(H), t, 2);
  hilbert_map(hecke_tau(H), hecke_t1t2(H), delta, prec);
  
  /* Sanity check: number of initialized elements */
  if (nb != hilbert_nb_cosets(ell))
    {
      flint_printf("(hecke_set_hilbert) Error: Hecke data structure initialized with %wd slots instead of the expected %wd\n", nb, hilbert_nb_cosets(ell));
      fflush(stdout);
      flint_abort();
    }
  
  if (v) flint_printf("(hecke_set_hilbert) Computing theta constants (%wd)", nb);
  fflush(stdout);

  /* Loop over all cosets to compute desired data */
  for (k = 0; k < nb; k++)
    {
      /* Print some status info? */
      if (v)
	{
	  if ((k+1) % 100 == 0)
	    {
	      flint_printf("\n(hecke_set_hilbert) (%wd/%wd)", k+1, nb);
	    }
	  flint_printf("."); fflush(stdout);
	}
      
      hilbert_coset(m, k, ell, delta);
      /* We know hilbert_coset consists only of scalars; in a next version,
	 change semantics of hilbert_cosets */
      fmpz_poly_mul(fmpz_poly_mat_entry(m, 1, 0),
		    fmpz_poly_mat_entry(m, 1, 0), beta);
      fmpz_poly_mul(fmpz_poly_mat_entry(m, 1, 1),
		    fmpz_poly_mat_entry(m, 1, 1), beta);
      hilbert_mat_map(gamma, m, delta);

      res = hecke_set_entry(H, k, gamma, prec);      
      if (!res)
	{
	  flint_printf("(hecke_set_hilbert) Warning: computation aborted due to low precision\n");
	  break;
	}
    }
  
  if (v) flint_printf("\n");

  fmpz_poly_mat_clear(m);
  fmpz_mat_clear(gamma);
  return res;  
}
