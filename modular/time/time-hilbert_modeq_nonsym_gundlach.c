
#include <stdio.h>
#include <flint/profiler.h> /* Flint profiler */
#include "modular.h"

/* Directory where data is written is TIMEDIR */

int main()
{
  slong ell;
  flint_rand_t state;
  FILE* data5;
  slong ell_max = 200;
  slong delta;
  slong k;
  fmpz_poly_t beta;
  
  fmpz_poly_struct num_vec[2];
  fmpz_t den;
  fmpq* rs;
  slong rs_bits = 4;
  timeit_t time;

  fmpz_poly_init(beta);
  for (k = 0; k < 2; k++) fmpz_poly_init(&num_vec[k]);
  fmpz_init(den);
  rs = _fmpq_vec_init(2);
  
  data5 = fopen(TIMEDIR "/data-hilbert_modeq_nonsym_gundlach", "w");
  
  flint_randinit(state);
  
  /* First line format: xmin xmax ymin ymax xlabel ylabel 
     or: xlabel ylabel */
  flint_fprintf(data5, "ell time(s)\n");

  delta = 5;
  for (k = 0; k < 2; k++)
    {
      fmpq_randbits(&rs[k], state, rs_bits);
    }
  
  for (ell = 2; ell < ell_max; ell++)
    {
      if (n_is_prime(ell) && hilbert_splits(beta, ell, delta))
	{
	  flint_printf("\n(time-hilbert_modeq_nonsym_gundlach) delta = %wd, ell = %wd, rs_bits = %wd\n", delta, ell, rs_bits);
	  timeit_start(time);      
	  hilbert_modeq_nonsym_gundlach_eval_Q(num_vec, den, rs, ell, beta, delta);
	  timeit_stop(time);
	  flint_fprintf(data5, "%wd %lf\n", ell, (double) time->cpu / 1000);
	}
    }

  fmpz_poly_clear(beta);
  for (k = 0; k < 2; k++) fmpz_poly_clear(&num_vec[k]);
  fmpz_clear(den);
  _fmpq_vec_clear(rs, 2);
  
  fclose(data5);
  flint_randclear(state);
  flint_cleanup();
  flint_printf("Done\n");
  return EXIT_SUCCESS;
}
