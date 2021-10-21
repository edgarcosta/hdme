
#include <stdio.h>
#include <profiler.h> /* Flint profiler */
#include "modular.h"

/* Directory where data is written is TIMEDIR */

int main()
{
  slong ell;
  flint_rand_t state;
  FILE* data5;
  FILE* data8;
  FILE* data13;
  FILE* data17;
  slong ell_max[4] = {1000, 300, 200, 100};
  slong delta;
  slong k;
  fmpz_poly_t beta;
  
  fmpz_poly_t num1, num2, num3;
  fmpz_t den;
  fmpq* rs;
  slong rs_bits = 10;
  timeit_t time;

  fmpz_poly_init(beta);
  fmpz_poly_init(num1);
  fmpz_poly_init(num2);
  fmpz_poly_init(num3);
  fmpz_init(den);
  rs = _fmpq_vec_init(2);
  
  data5 = fopen(TIMEDIR "/data-hilbert_modeq_eval_5", "w");
  data8 = fopen(TIMEDIR "/data-hilbert_modeq_eval_8", "w");
  data13 = fopen(TIMEDIR "/data-hilbert_modeq_eval_13", "w");
  data17 = fopen(TIMEDIR "/data-hilbert_modeq_eval_17", "w");
  
  flint_randinit(state);
  
  /* First line format: xmin xmax ymin ymax xlabel ylabel 
     or: xlabel ylabel */
  flint_fprintf(data5, "ell time(s)\n");
  flint_fprintf(data8, "ell time(s)\n");
  flint_fprintf(data13, "ell time(s)\n");
  flint_fprintf(data17, "ell time(s)\n");

  delta = 5;
  for (ell = 2; ell < ell_max[0]; ell++)
    {
      if (n_is_prime(ell) && hilbert_splits(beta, ell, delta))
	{
	  for (k = 0; k < 2; k++)
	    {
	      fmpq_randbits(&rs[k], state, rs_bits);
	    }
	  flint_printf("\n(time-hilbert_modeq_eval) delta = %wd, ell = %wd, rs_bits = %wd\n", delta, ell, rs_bits);
	  timeit_start(time);      
	  hilbert_modeq_sym_igusa_eval_Q(num1, num2, num3, den, rs, ell, delta);
	  timeit_stop(time);
	  flint_fprintf(data5, "%wd %lf\n", ell, (double) time->cpu / 1000);
	}
    }

  delta = 8;
  for (ell = 2; ell < ell_max[1]; ell++)
    {
      if (n_is_prime(ell) && hilbert_splits(beta, ell, delta))
	{
	  for (k = 0; k < 2; k++)
	    {
	      fmpq_randbits(&rs[k], state, rs_bits);
	    }
	  flint_printf("\n(time-hilbert_modeq_eval) delta = %wd, ell = %wd, rs_bits = %wd\n", delta, ell, rs_bits);
	  timeit_start(time);      
	  hilbert_modeq_sym_igusa_eval_Q(num1, num2, num3, den, rs, ell, delta);
	  timeit_stop(time);
	  flint_fprintf(data8, "%wd %lf\n", ell, (double) time->cpu / 1000);
	}
    }

  delta = 13;
  for (ell = 2; ell < ell_max[2]; ell++)
    {
      if (n_is_prime(ell) && hilbert_splits(beta, ell, delta))
	{
	  for (k = 0; k < 2; k++)
	    {
	      fmpq_randbits(&rs[k], state, rs_bits);
	    }
	  flint_printf("\n(time-hilbert_modeq_eval) delta = %wd, ell = %wd, rs_bits = %wd\n", delta, ell, rs_bits);
	  timeit_start(time);      
	  hilbert_modeq_sym_igusa_eval_Q(num1, num2, num3, den, rs, ell, delta);
	  timeit_stop(time);
	  flint_fprintf(data13, "%wd %lf\n", ell, (double) time->cpu / 1000);
	}
    }

  delta = 17;
  for (ell = 2; ell < ell_max[3]; ell++)
    {
      if (n_is_prime(ell) && hilbert_splits(beta, ell, delta))
	{
	  for (k = 0; k < 2; k++)
	    {
	      fmpq_randbits(&rs[k], state, rs_bits);
	    }
	  flint_printf("\n(time-hilbert_modeq_eval) delta = %wd, ell = %wd, rs_bits = %wd\n", delta, ell, rs_bits);
	  timeit_start(time);      
	  hilbert_modeq_sym_igusa_eval_Q(num1, num2, num3, den, rs, ell, delta);
	  timeit_stop(time);
	  flint_fprintf(data17, "%wd %lf\n", ell, (double) time->cpu / 1000);
	}
    }

  fmpz_poly_clear(beta);
  fmpz_poly_clear(num1);
  fmpz_poly_clear(num2);
  fmpz_poly_clear(num3);
  fmpz_clear(den);
  _fmpq_vec_clear(rs, 2);
  
  fclose(data5);
  fclose(data8);
  fclose(data13);
  fclose(data17);
  flint_randclear(state);
  flint_cleanup();
  flint_printf("Done\n");
  return EXIT_SUCCESS;
}
