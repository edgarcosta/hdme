
#include "modular.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("hilbert_modeq_gundlach_den....");
  fflush(stdout);

  flint_randinit(state);

  for (iter = 0; iter < 1 * arb_test_multiplier(); iter++)
    {
      acb_ptr t;
      fmpz_poly_mat_t m;
      slong ell = 11;
      fmpz_poly_t beta;
      fmpz_poly_t betabar;
      slong delta = 5;
      slong m_bits = 1 + n_randint(state, 5);
      slong prec = 500 + n_randint(state, 1000);
      slong n = hilbert_nb_cosets(ell, delta);
      int res;
      slong k;

      acb_ptr th2_vec;
      acb_ptr I_vec;
      acb_ptr stardets;
      acb_t scal;
      acb_t temp;
      acb_t den, den_test;
      acb_t star;
      
      t = _acb_vec_init(2);
      fmpz_poly_mat_init(m, 2, 2);
      fmpz_poly_init(beta);
      fmpz_poly_init(betabar);
      th2_vec = _acb_vec_init(32*n);
      I_vec = _acb_vec_init(8*n);
      stardets = _acb_vec_init(2*n);
      acb_init(scal);
      acb_init(temp);
      acb_init(den);
      acb_init(den_test);
      acb_init(star);

      hilbert_splits(beta, ell, delta);
      hilbert_conjugate(betabar, beta, delta);
      hilbert_halfspace_randtest(t, state, prec);
      hilbert_transform_randtest(m, state, m_bits);

      /* Denominator should be a modular form */
      res = hilbert_modeq_theta2_star(th2_vec, stardets, t, beta,
				      ell, delta, prec)
	&& hilbert_modeq_theta2_star(&th2_vec[16*n], &stardets[n], t, betabar,
				     ell, delta, prec);
      if (!res)
	{
	  flint_printf("FAIL (theta constants no.2)\n");
	  fflush(stdout);
	  flint_abort();
	}      
      modeq_cov(I_vec, th2_vec, n, prec);
      modeq_cov(&I_vec[4*n], &th2_vec[16*n], n, prec);

      acb_one(scal);
      for (k = 0; k < 2*n; k++)
	{
	  acb_pow_si(temp, &stardets[k], -10, prec);
	  acb_mul(scal, scal, temp, prec);
	}      
      hilbert_modeq_gundlach_den(den, I_vec, &I_vec[4*n], scal, ell, delta, prec);
            
      hilbert_star(star, m, t, delta, prec);
      acb_pow_ui(star, star, 2*10*n, prec);
      acb_mul(den, den, star, prec);
      
      /* Now compare with value at m*(t1,t2) */
      hilbert_transform(t, m, t, delta, prec);
      res = hilbert_modeq_theta2_star(th2_vec, stardets, t, beta,
				      ell, delta, prec)
	&& hilbert_modeq_theta2_star(&th2_vec[16*n], &stardets[n], t, betabar,
				     ell, delta, prec);
      if (!res)
	{
	  flint_printf("FAIL (theta constants no.2)\n");
	  fflush(stdout);
	  flint_abort();
	}      
      modeq_cov(I_vec, th2_vec, n, prec);
      modeq_cov(&I_vec[4*n], &th2_vec[16*n], n, prec);

      acb_one(scal);
      for (k = 0; k < 2*n; k++)
	{
	  acb_pow_si(temp, &stardets[k], -10, prec);
	  acb_mul(scal, scal, temp, prec);
	}      
      hilbert_modeq_gundlach_den(den_test, I_vec, &I_vec[4*n], scal, ell, delta, prec);
      
      if (!acb_overlaps(den, den_test))
	{
	  flint_printf("FAIL (overlap)\n");
	  acb_printd(den, 30); flint_printf("\n");
	  acb_printd(den_test, 30); flint_printf("\n");
	  fflush(stdout);
	  flint_abort();
	}

      _acb_vec_clear(t, 2);
      fmpz_poly_mat_clear(m);
      fmpz_poly_clear(beta);
      fmpz_poly_clear(betabar);
      _acb_vec_clear(th2_vec, 32*n);
      _acb_vec_clear(I_vec, 8*n);
      _acb_vec_clear(stardets, 2*n);
      acb_clear(scal);
      acb_clear(temp);
      acb_clear(den);
      acb_clear(den_test);
      acb_clear(star);
    }

  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

  
