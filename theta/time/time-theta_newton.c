
#include <stdio.h>
#include <flint/profiler.h> /* Flint profiler */
#include "theta.h"

#define TIME_THETA_NEWTON_ITER 100
#define TIME_THETA_NEWTON_MAXPREC 20000
/* Directory where data is written is TIMEDIR */

int main()
{
  slong iter;
  flint_rand_t state;
  FILE* data;

  data = fopen(TIMEDIR "/data-theta_newton", "w");
  flint_printf("time-theta_newton (%wd)", TIME_THETA_NEWTON_ITER);
  fflush(stdout);

  flint_randinit(state);

  /* First line format: xmin xmax ymin ymax xlabel ylabel 
     or: xlabel ylabel */
  flint_fprintf(data, "precision time(ms)\n");

  for (iter = 0; iter < TIME_THETA_NEWTON_ITER; iter++)
    {
      slong g = 2;
      slong prec = 10 + n_randint(state, TIME_THETA_NEWTON_MAXPREC);
      slong tol_bits = prec / 5;
      acb_mat_t tau;
      fmpz_mat_t m;
      acb_ptr th2;
      arb_t tol;
      
      timeit_t time;
      slong nb_iter;

      acb_mat_init(tau, g, g);
      fmpz_mat_init(m, 2*g, 2*g);
      th2 = _acb_vec_init(n_pow(2, 2*g));
      arb_init(tol);

      arb_set_si(tol, 1);
      arb_mul_2exp_si(tol, tol, -tol_bits); /* tol is larger than 2^(-prec) */

      siegel_halfspace_randtest(tau, state, prec);
      siegel_fundamental_domain(tau, m, tau, tol, prec);

      if (theta_use_newton(tau, prec)) /* tau lies in the right compact set */
	{
	  /* Do the computation */
	  TIMEIT_REPEAT(time, nb_iter);
	  theta2_newton(th2, tau, prec);
	  TIMEIT_END_REPEAT(time, nb_iter);
      
	  /* Write */
	  flint_fprintf(data, "%wd %lf\n", prec, (double) time->cpu / (double) nb_iter);
	  
	  flint_printf(".");
	  fflush(stdout);
	}
      else
	{
	  iter -= 1;
	}

      acb_mat_clear(tau);
      fmpz_mat_clear(m);
      _acb_vec_clear(th2, n_pow(2, 2*g));
      arb_clear(tol);
    }

  fclose(data);
  flint_randclear(state);
  flint_cleanup();
  flint_printf("done\n");
  return EXIT_SUCCESS;
}
