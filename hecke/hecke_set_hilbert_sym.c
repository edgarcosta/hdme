
#include "hilbert.h"

int hecke_set_hilbert_sym(hecke_t H, acb_srcptr t, const fmpz_poly_t beta,
			  slong ell, slong delta, slong prec)
{
  slong k;
  fmpz_poly_mat_t m;
  fmpz_mat_t gamma;
  fmpz_poly_t betabar;
  slong nb = hecke_nb(H);
  int res;
  int v = HECKE_VERBOSE;

  fmpz_poly_mat_init(m, 2, 2);
  fmpz_mat_init(gamma, 4, 4);
  fmpz_poly_init(betabar);

  /* Set Hecke context */
  hecke_ell(H) = ell;
  fmpz_poly_set(hecke_beta(H), beta);
  _acb_vec_set(hecke_t1t2(H), t, 2);
  hilbert_map(hecke_tau(H), hecke_t1t2(H), delta, prec);
  hilbert_conjugate(betabar, beta, delta);
  
  /* Sanity check: number of initialized elements */
  if (nb != 2*hilbert_nb_cosets(ell))
    {
      flint_printf("(hecke_set_hilbert_sym) Error: Hecke data structure initialized with %wd slots instead of the expected %wd\n", nb, 2*hilbert_nb_cosets(ell));
      fflush(stdout);
      flint_abort();
    }
  
  if (v) flint_printf("(hecke_set_hilbert_sym) Computing theta constants (%wd)", nb);
  fflush(stdout);
  
  /* Loop over all cosets to compute desired data */
  for (k = 0; k < nb; k++)
    {
      /* Print some status info? */
      if (v)
	{
	  if ((k+1) % 100 == 0)
	    {
	      flint_printf("\n(hecke_set_hilbert_sym) (%wd/%wd)", k+1, nb);
	    }
	  flint_printf("."); fflush(stdout);
	}

      if (k < hilbert_nb_cosets(ell))
	{
	  hilbert_coset(m, k, ell, delta);
	  /* We know hilbert_coset consists only of scalars; in a next version,
	     change semantics of hilbert_cosets */
	  fmpz_poly_mul(fmpz_poly_mat_entry(m, 1, 0),
			fmpz_poly_mat_entry(m, 1, 0), beta);
	  fmpz_poly_mul(fmpz_poly_mat_entry(m, 1, 1),
			fmpz_poly_mat_entry(m, 1, 1), beta);
	  hilbert_mat_map(gamma, m, delta);
	}
      else
	{
	  hilbert_coset(m, k - hilbert_nb_cosets(ell), ell, delta);
	  fmpz_poly_mul(fmpz_poly_mat_entry(m, 1, 0),
			fmpz_poly_mat_entry(m, 1, 0), betabar);
	  fmpz_poly_mul(fmpz_poly_mat_entry(m, 1, 1),
			fmpz_poly_mat_entry(m, 1, 1), betabar);
	  hilbert_mat_map(gamma, m, delta);	  
	}

      res = hecke_set_entry(H, k, gamma, prec);      
      if (!res)
	{
	  flint_printf("(hecke_set_hilbert_sym) Warning: computation aborted due to low precision\n");
	  break;
	}
    }
  
  if (v) flint_printf("\n");

  fmpz_poly_mat_clear(m);
  fmpz_mat_clear(gamma);
  fmpz_poly_clear(betabar);
  return res;  
}

