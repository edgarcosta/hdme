
#include "hecke.h"

/* Test memory allocation only, not contents */

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("hecke_set....");
  fflush(stdout);

  flint_randinit(state);
  
  for (iter = 0; iter < 10 * arb_test_multiplier(); iter++)
    {
      acb_mat_t tau;
      fmpz* I;
      acb_ptr t;
      hecke_t H;
      slong delta;
      slong ell;
      fmpz_poly_t beta;
      slong prec = 500 + n_randint(state, 1000);

      acb_mat_init(tau, 2, 2);
      I = _fmpz_vec_init(4);
      t = _acb_vec_init(2);
      fmpz_poly_init(beta);

      /* Test setting context data */
      hecke_init(H, 1);
      
      fmpz_set(&I[0], 1+n_randint(state, 100));
      fmpz_set(&I[1], 1+n_randint(state, 100));
      fmpz_set(&I[2], 1+n_randint(state, 100));
      fmpz_set(&I[3], 1+n_randint(state, 100));
      hecke_set_I(H, I, prec);

      siegel_fundamental_domain_randtest(tau, state, prec);
      hecke_set_tau(H, tau, prec);

      hilbert_halfspace_randtest(t, state, prec);
      hecke_set_t1t2(H, t, prec);
      
      hecke_clear(H);

      /* Test setting Hecke correspondences */
      ell = 3;
      hecke_init(H, siegel_nb_cosets(ell));
      hecke_set_siegel(H, tau, ell, prec);
      hecke_clear(H);

      delta = 5;
      ell = 11;      
      hilbert_splits(beta, ell, delta);
      hecke_init(H, hilbert_nb_cosets(ell));
      hecke_set_hilbert(H, t, beta, ell, delta, prec);
      hecke_clear(H);

      hecke_init(H, 2*hilbert_nb_cosets(ell));
      hecke_set_hilbert_sym(H, t, beta, ell, delta, prec);
      hecke_clear(H);

      ell = 2;
      hecke_init(H, siegel_nb_T1_cosets(ell));
      hecke_set_t1(H, tau, ell, prec);
      hecke_clear(H);
    }

  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}
