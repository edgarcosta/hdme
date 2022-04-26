
#include <stdio.h>
#include <math.h>
#include <flint/profiler.h> /* Flint profiler */
#include "modular.h"

#define TIME_HILBERT_NB_PRIMES 100
#define TIME_HILBERT_BITS 8
/* Directory where data is written is TIMEDIR */

int main()
{
  slong iter;
  flint_rand_t state;
  FILE* data;
  slong j_bits = TIME_HILBERT_BITS;  
  fmpq* j;
  fmpz_t j_num, j_den;
  slong delta = 5;
  slong k;
  
  flint_randinit(state);
  
  j = _fmpq_vec_init(2);
  fmpz_init(j_num);
  fmpz_init(j_den);

  /* Change starting values */
  fmpz_randbits(j_num, state, 100);
  
  for (k = 0; k < 2; k++)
    {
      fmpz_randbits(j_num, state, j_bits);
      fmpz_randbits(j_den, state, j_bits);
      /*fmpz_one(j_den);*/
      fmpq_set_fmpz_frac(&j[k], j_num, j_den);
    }
  
  flint_printf("(time-hilbert_gundlach) Evaluate Hilbert modular equations at Gundlach invariants\n");
  for (k = 0; k < 2; k++)
    {
      fmpq_print(&j[k]); flint_printf("\n");
    }

  /* First line format: xmin xmax ymin ymax xlabel ylabel 
     or: xlabel ylabel */
  data = fopen(TIMEDIR "/data-hilbert_gundlach", "w");
  flint_fprintf(data, "0 %wd 0 200 l time(s)\n", n_nth_prime(TIME_HILBERT_NB_PRIMES)+1);
  fclose(data);
  
  for (iter = 1; iter <= TIME_HILBERT_NB_PRIMES; iter++)
    {
      fmpz_poly_struct num_vec[2];
      fmpz_t den;
      slong ell = n_nth_prime(iter);
      timeit_t time;
      int res;
      
      for (k = 0; k < 2; k++) fmpz_poly_init(&num_vec[k]);
      fmpz_init(den);

      flint_printf("\n(time-hilbert_gundlach) ell = %wd\n", ell);
      timeit_start(time);      
      res = hilbert_modeq_gundlach_eval_Q(num_vec, den, j, ell, delta);
      timeit_stop(time);
      if (res == 1)
	{
	  data = fopen(TIMEDIR "/data-hilbert_gundlach", "a");
	  flint_fprintf(data, "%wd %lf\n", ell,
			(double) time->cpu / 1000);
	  fclose(data);
	}

      for (k = 0; k < 2; k++) fmpz_poly_clear(&num_vec[k]);
      fmpz_clear(den);
    }
  
  _fmpq_vec_clear(j, 2);
  fmpz_clear(j_num);
  fmpz_clear(j_den);

  flint_randclear(state);
  flint_cleanup();
  flint_printf("Done\n");
  return EXIT_SUCCESS;
}
