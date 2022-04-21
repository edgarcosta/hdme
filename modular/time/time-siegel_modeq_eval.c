
#include <stdio.h>
#include <flint/profiler.h> /* Flint profiler */
#include "modular.h"

#define TIME_SIEGEL_MODEQ_EVAL_ITER 4
#define TIME_SIEGEL_MODEQ_EVAL_STEP 5 /* Up to 25 bits in length */ 
/* Directory where data is written is TIMEDIR */

int main()
{
  slong iter;
  flint_rand_t state;
  FILE* data2;
  FILE* data3;
  FILE* data5;
  FILE* data7;
  FILE* data11;
  slong j_bits = 0;
  
  data2 = fopen(TIMEDIR "/data-siegel_modeq_eval_2", "w");
  data3 = fopen(TIMEDIR "/data-siegel_modeq_eval_3", "w");
  data5 = fopen(TIMEDIR "/data-siegel_modeq_eval_5", "w");
  data7 = fopen(TIMEDIR "/data-siegel_modeq_eval_7", "w");
  data11 = fopen(TIMEDIR "/data-siegel_modeq_eval_11", "w");
  flint_printf("time-siegel_modeq_eval (%wd)\n", TIME_SIEGEL_MODEQ_EVAL_ITER);
  fflush(stdout);

  flint_randinit(state);

  /* First line format: xmin xmax ymin ymax xlabel ylabel 
     or: xlabel ylabel */
  flint_fprintf(data2, "bitlength time(s)\n");
  flint_fprintf(data3, "bitlength time(s)\n");
  flint_fprintf(data5, "bitlength time(s)\n");
  flint_fprintf(data7, "bitlength time(s)\n");
  flint_fprintf(data11, "bitlength time(s)\n");

  for (iter = 0; iter < TIME_SIEGEL_MODEQ_EVAL_ITER; iter++)
    {
      fmpz_poly_struct num_vec[3];
      fmpz_t den;
      fmpq* j;
      fmpz_t j_num, j_den;
      slong ell;
      slong k;
      timeit_t time;
      
      j_bits += TIME_SIEGEL_MODEQ_EVAL_STEP;

      for (k = 0; k < 3; k++) fmpz_poly_init(&num_vec[k]);
      fmpz_init(den);
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

      ell = 2;
      flint_printf("\n(time-siegel_modeq_eval) ell = %wd, j_bits = %wd\n", ell, j_bits);
      timeit_start(time);      
      siegel_modeq_eval_Q(num_vec, den, j, ell);
      timeit_stop(time);
      flint_fprintf(data2, "%wd %lf\n", j_bits, (double) time->cpu / 1000);

      ell = 3;
      flint_printf("\n(time-siegel_modeq_eval) ell = %wd, j_bits = %wd\n", ell, j_bits);
      timeit_start(time);      
      siegel_modeq_eval_Q(num_vec, den, j, ell);
      timeit_stop(time);
      flint_fprintf(data3, "%wd %lf\n", j_bits, (double) time->cpu / 1000);

      /*
      ell = 5;
      flint_printf("\n(time-siegel_modeq_eval) ell = %wd, j_bits = %wd\n", ell, j_bits);
      timeit_start(time);      
      siegel_modeq_eval_Q(num_vec, den, j, ell);
      timeit_stop(time);
      flint_fprintf(data5, "%wd %lf\n", j_bits, (double) time->cpu / 1000);
      
      
      ell = 7;
      flint_printf("\n(time-siegel_modeq_eval) ell = %wd, j_bits = %wd\n", ell, j_bits);
      timeit_start(time);      
      siegel_modeq_eval_Q(num_vec, den, j, ell);
      timeit_stop(time);
      flint_fprintf(data7, "%wd %lf\n", prec, (double) time->cpu / 1000);
      
      ell = 11;
      flint_printf("\n(time-siegel_modeq_eval) ell = %wd, j_bits = %wd\n", ell, j_bits);
      timeit_start(time);      
      siegel_modeq_eval_Q(num_vec, den, j, ell);
      timeit_stop(time);
      flint_fprintf(data11, "%wd %lf\n", prec, (double) time->cpu / 1000);
      */

      for (k = 0; k < 3; k++) fmpz_poly_clear(&num_vec[k]);
      fmpz_clear(den);
      _fmpq_vec_clear(j, 3);
      fmpz_clear(j_num);
      fmpz_clear(j_den);
    }

  fclose(data2);
  fclose(data3);
  fclose(data5);
  fclose(data7);
  fclose(data11);
  flint_randclear(state);
  flint_cleanup();
  flint_printf("Done\n");
  return EXIT_SUCCESS;
}
