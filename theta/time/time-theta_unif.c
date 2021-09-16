
#include <stdio.h>
#include <profiler.h> /* Flint profiler */
#include "theta.h"

#define TIME_THETA_UNIF_ITER 200
#define TIME_THETA_UNIF_MAXPREC 50000
#define TIME_THETA_UNIF_MAXMULT 20
/* Directory where data is written is TIMEDIR */

int main()
{
  slong iter;
  flint_rand_t state;
  FILE* data;
  timeit_t total_time;

  data = fopen(TIMEDIR "/data-theta_unif", "w");
  flint_printf("time-theta_unif (%wd)", TIME_THETA_UNIF_ITER);
  fflush(stdout);

  flint_randinit(state);
  timeit_start(total_time);

  /* First line format: xmin xmax ymin ymax xlabel ylabel 
     or: xlabel ylabel */
  flint_fprintf(data, "precision time(ms)\n");

  for (iter = 0; iter < TIME_THETA_UNIF_ITER; iter++)
    {
      slong g = 2;
      slong prec = 10 + n_randint(state, TIME_THETA_UNIF_MAXPREC);
      slong tol_bits = prec / 5;
      slong mult = 1 + n_randint(state, TIME_THETA_UNIF_MAXMULT);
      acb_mat_t tau;
      sp2gz_t m;
      acb_ptr th2;
      arb_t tol;
      
      timeit_t time;
      slong nb_iter;

      acb_mat_init(tau, g, g);
      sp2gz_init(m, g);
      th2 = _acb_vec_init(n_pow(2, 2*g));
      arb_init(tol);

      arb_set_si(tol, 1);
      arb_mul_2exp_si(tol, tol, -tol_bits); /* tol is larger than 2^(-prec) */

      siegel_halfspace_randtest(tau, state, prec);
      siegel_fundamental_domain(tau, m, tau, tol, prec);
      acb_mat_scalar_mul_si(tau, tau, mult, prec);
      siegel_fundamental_domain(tau, m, tau, tol, prec);

      
      /* Do the computation */
      TIMEIT_REPEAT(time, nb_iter);
      theta2_unif(th2, tau, prec);
      TIMEIT_END_REPEAT(time, nb_iter);
      
      /* Write */
      flint_fprintf(data, "%wd %lf\n", prec, (double) time->cpu / (double) nb_iter);
      
      flint_printf(".");
      fflush(stdout);
  
      acb_mat_clear(tau);
      sp2gz_clear(m);
      _acb_vec_clear(th2, n_pow(2, 2*g));
      arb_clear(tol);
    }

  timeit_stop(total_time);
  fclose(data);
  flint_randclear(state);
  flint_cleanup();
  flint_printf("done. Total time: %wdms\n", total_time->cpu);
  return EXIT_SUCCESS;
}
