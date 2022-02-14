
#include <stdio.h>
#include <math.h>
#include <profiler.h> /* Flint profiler */
#include "modular.h"

#define TIME_SIEGEL_NB_PRIMES 3 /* 6 */
#define TIME_SIEGEL_BITS 8
/* Directory where data is written is TIMEDIR */

int main()
{
  slong iter;
  flint_rand_t state;
  FILE* data;
  slong j_bits = TIME_SIEGEL_BITS;  
  fmpq* j;
  fmpz_t j_num, j_den;
  slong k;
  
  data = fopen(TIMEDIR "/data-siegel", "w");
  flint_randinit(state);
  
  j = _fmpq_vec_init(3);
  fmpz_init(j_num);
  fmpz_init(j_den);
  
  for (k = 0; k < 3; k++)
    {
      fmpz_randbits(j_num, state, j_bits);
      fmpz_randbits(j_den, state, j_bits);
      /*fmpz_one(j_den);*/
      fmpq_set_fmpz_frac(&j[k], j_num, j_den);
    }
  
  flint_printf("(time-siegel) Evaluate Siegel modular equations at Igusa invariants\n");
  for (k = 0; k < 3; k++)
    {
      fmpq_print(&j[k]); flint_printf("\n");
    }

  /* First line format: xmin xmax ymin ymax xlabel ylabel 
     or: xlabel ylabel */
  flint_fprintf(data, "0 %wd 0 20 l time(s)\n", n_nth_prime(TIME_SIEGEL_NB_PRIMES)+1);
  
  for (iter = 1; iter <= TIME_SIEGEL_NB_PRIMES; iter++)
    {
      fmpz_poly_struct num_vec[3];
      fmpz_t den;
      slong ell = n_nth_prime(iter);
      timeit_t time;
      
      for (k = 0; k < 3; k++) fmpz_poly_init(&num_vec[k]);
      fmpz_init(den);

      flint_printf("\n(time-siegel) ell = %wd\n", ell);
      timeit_start(time);      
      siegel_modeq_eval_Q(num_vec, den, j, ell);
      timeit_stop(time);
      flint_fprintf(data, "%wd %lf\n", ell,
		    (double) time->cpu / 1000);

      for (k = 0; k < 3; k++) fmpz_poly_clear(&num_vec[k]);
      fmpz_clear(den);
    }
  
  _fmpq_vec_clear(j, 3);
  fmpz_clear(j_num);
  fmpz_clear(j_den);

  fclose(data);
  flint_randclear(state);
  flint_cleanup();
  flint_printf("Done\n");
  return EXIT_SUCCESS;
}
