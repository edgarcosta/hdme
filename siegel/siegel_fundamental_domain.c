
#include "siegel.h"

int siegel_fundamental_domain(acb_mat_t w, sp2gz_t m,
			      const acb_mat_t z, const arb_t tol, slong prec)
{
  int stop = 0;
  int res = 1;
  int j;
  
  sp2gz_t test_matrix;
  sp2gz_t step;
  fmpz_mat_t u;
  arb_mat_t im;
  acb_mat_t x;
  acb_mat_t star;
  acb_t det;
  arb_t absdet;

  sp2gz_init(test_matrix, m->g);
  sp2gz_init(step, m->g);
  fmpz_mat_init(u, m->g, m->g);
  arb_mat_init(im, m->g, m->g);
  acb_mat_init(x, m->g, m->g);
  acb_mat_init(star, m->g, m->g);
  acb_init(det);
  arb_init(absdet);
  
  sp2gz_one(m);

  while (!stop)
    {
      /* Update x */
      res = res && siegel_transform(x, m, z, prec);
      if (!res)
	{
	  flint_printf("\n(siegel_fundamental_domain) Precision lost in update:\n");
	  flint_printf("m = "); sp2gz_print(m); flint_printf("\n\n");
	  flint_printf("z = "); acb_mat_printd(z, 30); flint_printf("\n\n");
	  flint_printf("tol = "); arb_printd(tol, 5); flint_printf("\n");
	  break;
	}
      
      /* Reduce imaginary part */
      acb_mat_get_imag(im, x);
      res = arb_mat_minkowski_reduce(im, u, im, tol, prec);
      sp2gz_set_diagonal(step, u);
      
      sp2gz_mul(m, step, m);
      res = res && siegel_transform(x, m, z, prec);
      if (!res)
	{
	  flint_printf("\n(siegel_fundamental_domain) Precision lost in Minkowski reduction:\n");
	  flint_printf("im = "); arb_mat_printd(im, 30); flint_printf("\n\n");
	  flint_printf("u = "); fmpz_mat_print_pretty(u); flint_printf("\n\n");
	  flint_printf("tol = "); arb_printd(tol, 5); flint_printf("\n");
	  break;
	}
      
      /* Reduce real part */
      res = siegel_reduce_real(x, step, x, tol, prec);
      sp2gz_mul(m, step, m);
      if (!res)
	{
	  flint_printf("\n(siegel_fundamental_domain) Precision lost in real reduction:\n");
	  flint_printf("x = "); acb_mat_printd(x, 30); flint_printf("\n\n");
	  flint_printf("tol = "); arb_printd(tol, 5); flint_printf("\n");
	  break;
	}

      /* Apply test matrix */
      /* To be improved: compute determinants at smaller precision,
	 select the smallest one */
      stop = 1;
      for (j = 0; j < siegel_nb_test_matrices(m->g); j++)
	{
	  siegel_test_matrix(step, j);
	  siegel_star(star, step, x, prec);
	  acb_mat_det(det, star, prec);
	  acb_abs(absdet, det, prec);
	  arb_sub_si(absdet, absdet, 1, prec);
	  if (arb_is_negative(absdet))
	    {
	      stop = 0;
	      sp2gz_mul(m, step, m);
	      break; /* for loop */
	    }
	}
    }
  
  res = res && siegel_transform(w, m, z, prec);

  sp2gz_clear(test_matrix);
  sp2gz_clear(step);
  fmpz_mat_clear(u);
  arb_mat_clear(im);
  acb_mat_clear(x);
  acb_mat_clear(star);
  acb_clear(det);
  arb_clear(absdet);
  return res;
}
