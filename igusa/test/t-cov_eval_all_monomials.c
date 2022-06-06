
#include "igusa.h"

int main()
{
  slong iter;
  flint_rand_t state;
  
  flint_printf("cov_eval_all_monomials....");
  fflush(stdout);

  flint_randinit(state);

  for (iter = 0; iter < 100 * arb_test_multiplier(); iter++)
    {
      slong nb = 1 + n_randint(state, 10);
      slong mag_bits = 10;
      slong prec = 200;
      slong* weights;
      slong wt;
      acb_ptr I;
      acb_ptr ev;
      slong m;
      slong k;

      weights = flint_malloc(nb * sizeof(slong));
      I = _acb_vec_init(nb);      
      for (k = 0; k < nb; k++)
	{
	  weights[k] = 1 + n_randint(state, 10);
	  acb_randtest_precise(&I[k], state, prec, mag_bits);
	}
      wt = n_randint(state, 20);
      m = cov_nb_monomials(wt, nb, weights);
      ev = _acb_vec_init(m);

      _acb_vec_clear(ev, m);
      _acb_vec_clear(I, nb);
      flint_free(weights);			   
    }

  flint_randclear(state);
  flint_cleanup();
  flint_printf("PASS\n");
  return EXIT_SUCCESS;
}

      
