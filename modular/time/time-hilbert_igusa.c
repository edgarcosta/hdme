
#include <stdio.h>
#include <profiler.h> /* Flint profiler */
#include "modular.h"

/* Directory where data is written is TIMEDIR */

int main()
{
  slong ell;
  flint_rand_t state;
  FILE* data;
  slong ell_max = 50;
  slong delta;
  slong discriminants[30] = {5, 8, 12, 13, 17, 21, 24, 28, 29, 33, 37, 40, 41, 44, 53, 56, 57, 60, 61, 65, 69, 73, 76, 77, 85, 88, 89, 92, 93, 97};
  char filename[200];
  slong j, k;
  fmpz_poly_t beta;
  
  fmpz_poly_struct num_vec[3];
  fmpz_t den;
  fmpq* rs;
  slong rs_bits = 8;
  timeit_t time;

  fmpz_poly_init(beta);
  for (k = 0; k < 3; k++) fmpz_poly_init(&num_vec[k]);
  fmpz_init(den);
  rs = _fmpq_vec_init(2);
  
  flint_randinit(state);

  /* Change values */
  fmpq_randbits(&rs[0], state, 200);
  for (k = 0; k < 2; k++)
    {
      fmpq_randbits(&rs[k], state, rs_bits);
    }

  for (j = 0; j < 30; j++)
    {
      delta = discriminants[j];
      flint_sprintf(filename, TIMEDIR "/data-hilbert_igusa_%wd", delta);
      
      /* First line format: xmin xmax ymin ymax xlabel ylabel 
	 or: xlabel ylabel */
      data = fopen(filename, "w");
      flint_fprintf(data, "ell time(s)\n");
      fclose(data);
      
      flint_printf("(time-hilbert_igusa) Evaluate Hilbert modular equations for discriminant %wd at parameters\n", delta);
      for (k = 0; k < 2; k++)
	{
	  fmpq_print(&rs[k]); flint_printf("\n");
	}
      for (ell = 2; ell < ell_max; ell++)
	{
	  if (n_is_prime(ell) && hilbert_splits(beta, ell, delta))
	    {
	      flint_printf("\n(time-hilbert_igusa) ell = %wd\n", ell);
	      timeit_start(time);      
	      hilbert_modeq_igusa_eval_Q(num_vec, den, rs, ell, delta);
	      timeit_stop(time);
	      data = fopen(filename, "a");
	      flint_fprintf(data, "%wd %lf\n", ell, (double) time->cpu / 1000);
	      fclose(data);
	    }
	}
    }

  fmpz_poly_clear(beta);
  for (k = 0; k < 3; k++) fmpz_poly_clear(&num_vec[k]);
  fmpz_clear(den);
  _fmpq_vec_clear(rs, 2);
  
  flint_randclear(state);
  flint_cleanup();
  flint_printf("Done\n");
  return EXIT_SUCCESS;
}
