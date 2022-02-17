
#include <assert.h>
#include "flint.h"
#include "arb.h"
#include "modular.h"

int main()
{
  fmpz* I;
  fmpq* j;
  slong ell = 3;
  fmpz_poly_struct num_vec[3];
  fmpz_t den;
  fmpz_poly_struct factors[40]; /* At most degree of modular eqs of level 3 */
  slong exps[40];
  slong nb_factors = 0;
  fmpz_poly_t fac;
  acb_poly_t fac_acb;
  acb_ptr roots;
  acb_ptr igusa_tuples;
  acb_poly_struct nums_acb[12];
  acb_ptr dens_acb;
  acb_poly_t trace_num;
  acb_poly_t norm_num;
  fmpz_poly_t norm_rounded;
  acb_t rescale;
  arf_t max_radius;
  fmpq* ratl_roots;
  
  slong prec;
  slong k, i;
  int res;

  /* Init everything */
  I = _fmpz_vec_init(4);
  j = _fmpq_vec_init(3);
  for (k = 0; k < 3; k++) fmpz_poly_init(&num_vec[k]);
  fmpz_init(den);
  for (k = 0; k < 40; k++) fmpz_poly_init(&factors[k]);
  fmpz_poly_init(fac);
  acb_poly_init(fac_acb);
  roots = _acb_vec_init(4);
  igusa_tuples = _acb_vec_init(12);
  for (k = 0; k < 12; k++) acb_poly_init(&nums_acb[k]);
  dens_acb = _acb_vec_init(4);
  acb_poly_init(trace_num);
  acb_poly_init(norm_num);
  fmpz_poly_init(norm_rounded);
  acb_init(rescale);
  arf_init(max_radius);
  ratl_roots = _fmpq_vec_init(40);
  
  /* Set values of I according to LMFDB web page
     https://olive.lmfdb.xyz/Genus2Curve/Q/277/a/277/1
     https://beta.lmfdb.org/Genus2Curve/Q/277/a/277/1 */
  fmpz_set_si(&I[0], 64);
  fmpz_set_si(&I[1], 352);
  fmpz_set_si(&I[2], 9552);
  fmpz_set_si(&I[3], -1108);

  /* We need I6', not I6 */
  igusa_I6prime_fmpz(&I[2], I);
  igusa_from_cov_fmpz(j, I);
  
  /* Evaluate modular equations over QQ */
  siegel_modeq_eval_Q(num_vec, den, j, ell);
  modeq_factor_Q(&nb_factors, factors, exps, &num_vec[0]);
  
  /* Isolate factor of degree 4 */
  for (k = 0; k < nb_factors; k++)
    {
      if (fmpz_poly_degree(&factors[k]) == 4)
	{
	  fmpz_poly_set(fac, &factors[k]);
	  fmpz_poly_primitive_part(fac, fac);
	  flint_printf("Degree 4 factor:\n");
	  fmpz_poly_print_pretty(fac, "x");
	  flint_printf("\n");
	  break;
	}
    }

  /* Compute vector of four Igusa tuples (j1, j2, j3) where j1 is a
     root of factor */
  prec = 14000;
  acb_poly_set_fmpz_poly(fac_acb, fac, prec);
  nb_factors = acb_poly_find_roots(roots, fac_acb, NULL, 0, prec);
  assert (nb_factors == 4);
  flint_printf("Complex roots:\n");
  for (k = 0; k < 4; k++)
    {
      acb_printd(&roots[k], 30); flint_printf("\n");
    }

  
  flint_printf("\n(3,3)-isogenous Igusa invariants in number field of degree 4 embedded in CC:\n");
  for (k = 0; k < 4; k++)
    {
      res = modeq_isog_invariants_C(&igusa_tuples[3*k],
				    num_vec, &roots[k], 3, prec);
      flint_printf("[");
      acb_printd(&igusa_tuples[3*k], 10); flint_printf(",\n ");
      acb_printd(&igusa_tuples[3*k+1], 10); flint_printf(",\n ");
      acb_printd(&igusa_tuples[3*k+2], 10); flint_printf("]\n");
    }

  /* Sanity check: what's the trace of j1, j2, j3? */
  /* We find 277^{-6}*269^{-4}*(approx. of integer) each time. This
     can be proved directly, since j1, j2, j3 are known algebraically */

  /* Evaluate modular equations over C at all these four points */
  for (k = 0; k < 4; k++)
    {
      siegel_modeq_eval_C(&nums_acb[3*k], &dens_acb[k],
			  &igusa_tuples[3*k], 3, prec);
    }
  
  /* Set rescale factor: denominator for j1, j2, j3 as algebraic
     numbers */
  acb_set_si(rescale, 277*269);
  acb_pow_ui(rescale, rescale, 4, prec);
  /* Now take an appropriate power (degree of modular equations) */
  acb_pow_ui(rescale, rescale, 10*40/3, prec);      

  /* Compute trace of rescaled modular equations */
  for (k = 0; k < 3; k++)
    {
      acb_poly_zero(trace_num);
      for (i = 0; i < 4; i++)
	{
	  acb_poly_add(trace_num, trace_num, &nums_acb[k+3*i], prec);
	}
      acb_poly_scalar_mul(trace_num, trace_num, rescale, prec);
      res = modeq_round_poly(&num_vec[k], max_radius, trace_num, 40);
      if (!res)
      	{
      	  flint_printf("Could not recognize integer coefficients: max radius ");
      	  arf_printd(max_radius, 10);
      	  flint_printf("\nConstant coefficient:\n");
      	  arb_printd(acb_realref(acb_poly_get_coeff_ptr(trace_num, 0)), 5000);
      	  flint_printf("\n");
      	  flint_abort();
      	}
    }
  /* We can also multiply all conjugates, but need higher precision (around 50000) */
  /* acb_poly_one(norm_num);
  for (i = 0; i < 4; i++)
    {
      acb_poly_mul(norm_num, norm_num, &nums_acb[3*i], prec);
    }
  acb_pow_ui(rescale, rescale, 4, prec);
  acb_poly_scalar_mul(norm_num, norm_num, rescale, prec);
  res = modeq_round_poly(norm_rounded, max_radius, norm_num, 4*40);
  if (!res)
    {
      flint_printf("Could not recognize integer coefficients: max radius ");
      arf_printd(max_radius, 10);
      flint_printf("\n");
      flint_abort();
      } */
  
  /* Now num1 should have (at least) two rational roots: the j1 we
     started with, and j1 of the (9,3,3)-isogenous p.p.a.s. */
  modeq_roots_Q(&nb_factors, ratl_roots, exps, &num_vec[0]);
  assert (nb_factors >= 2);

  flint_printf(" We find two rational roots:\n");
  for (k = 0; k < nb_factors; k++)
    {
      fmpq_print(&ratl_roots[k]); flint_printf("\n");
    }
  /* A priori we just have a root of the trace; to check it's indeed a
     root, we could just reconstruct the coefficients of &num_vec[0]
     algebraically in the number field, and check that evaluation
     vanishes as it should. Or, take the norm. */
  flint_printf("Igusa invariants of first genus 2 curve were:\n");
  for (k = 0; k < 3; k++)
    {
      fmpq_print(&j[k]); flint_printf("\n");
    }
  /* Take other root */
  fmpq_set_si(&ratl_roots[0], -274862096729792, 277);
  modeq_isog_invariants_Q(j, num_vec, &ratl_roots[0], 3);
  flint_printf("Igusa invariants of (9,3,3)-isogenous genus 2 curve:\n");
  for (k = 0; k < 3; k++)
    {
      fmpq_print(&j[k]); flint_printf("\n");
    }

  flint_printf("This matches with the other LMFDB curve:\n");
  /* Web page is https://beta.lmfdb.org/Genus2Curve/Q/277/a/277/2 */
  fmpz_set_si(&I[0], 4480);
  fmpz_set_si(&I[1], 1370512);
  fmpz_set_si(&I[2], 1511819744);
  fmpz_set_si(&I[3], -1108);
  /* We need I6', not I6 */
  igusa_I6prime_fmpz(&I[2], I);
  igusa_from_cov_fmpz(j, I);
  for (k = 0; k < 3; k++)
    {
      fmpq_print(&j[k]); flint_printf("\n");
    }  
  
  /* Clear everything */  
  _fmpz_vec_clear(I, 4);
  _fmpq_vec_clear(j, 3);
  for (k = 0; k < 3; k++) fmpz_poly_clear(&num_vec[k]);
  fmpz_clear(den);
  for (k = 0; k < 40; k++) fmpz_poly_clear(&factors[k]);
  fmpz_poly_clear(fac);
  acb_poly_clear(fac_acb);
  _acb_vec_clear(roots, 4);
  _acb_vec_clear(igusa_tuples, 12);
  for (k =0; k < 12; k++) acb_poly_clear(&nums_acb[k]);
  _acb_vec_clear(dens_acb, 4);
  acb_poly_clear(trace_num);
  acb_poly_clear(norm_num);
  fmpz_poly_clear(norm_rounded);
  acb_clear(rescale);
  arf_clear(max_radius);
  _fmpq_vec_clear(ratl_roots, 40);

  flint_cleanup();
  return EXIT_SUCCESS;
}
