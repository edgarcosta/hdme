
#include <stdio.h>
#include <profiler.h> /* Flint profiler */
#include "theta.h"

#define TIME_THETA_NAIVE_ITER 50
#define TIME_THETA_NAIVE_MAXPREC 10000
/* Directory where data is written is TIMEDIR */

int main()
{
  slong iter;
  flint_rand_t state;
  FILE* data;

  data = fopen(TIMEDIR "/data-theta_naive", "w");
  flint_printf("time-theta_naive (%wd)", TIME_THETA_NAIVE_ITER);
  fflush(stdout);

  flint_randinit(state);

  /* First line format: xmin xmax ymin ymax xlabel ylabel 
     or: xlabel ylabel */
  flint_fprintf(data, "precision time(ms)\n");

  for (iter = 0; iter < TIME_THETA_NAIVE_ITER; iter++)
    {
      slong g = 2;
      slong prec = 10 + n_randint(state, TIME_THETA_NAIVE_MAXPREC);
      slong mult = 1; /*+ n_randint(state, 5);*/ /* Closer to the cusp */
      acb_mat_t tau;
      acb_ptr th2;
      arb_mat_t im;
      arb_t lambda;
      
      timeit_t time;
      slong nb_iter;

      acb_mat_init(tau, g, g);
      arb_mat_init(im, g, g);
      th2 = _acb_vec_init(n_pow(2, 2*g));
      arb_init(lambda);

      siegel_halfspace_randtest(tau, state, prec);
      acb_mat_scalar_mul_si(tau, tau, mult, prec);
      acb_mat_get_imag(im, tau);
      arb_mat_lambda(lambda, im, prec);
      /* Scale tau such that im(tau) has an eigenvalue sqrt(3)/4 */
      acb_mat_scalar_div_arb(tau, tau, lambda, prec);
      arb_sqrt_ui(lambda, 3, prec);
      arb_div_si(lambda, lambda, 4, prec);
      acb_mat_scalar_mul_arb(tau, tau, lambda, prec);

      /* Do the computation */
      TIMEIT_REPEAT(time, nb_iter);
      theta2_naive(th2, tau, prec);
      TIMEIT_END_REPEAT(time, nb_iter);

      /* Write */
      flint_fprintf(data, "%wd %lf\n", prec, (double) time->cpu / (double) nb_iter);
      
      acb_mat_clear(tau);
      arb_mat_clear(im);
      _acb_vec_clear(th2, n_pow(2, 2*g));
      arb_clear(lambda);

      flint_printf(".");
      fflush(stdout);
    }

  fclose(data);
  flint_randclear(state);
  flint_cleanup();
  flint_printf("done\n");
  return EXIT_SUCCESS;
}
