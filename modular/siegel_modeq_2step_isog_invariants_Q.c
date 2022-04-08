
#include <assert.h>
#include "modular.h"

/* Compute the list of primitive, irreducible polynomials dividing
   the nonzero polynomial poly, of degree dividing l+1, and not linear */
static void
siegel_modeq_2step_acceptable_factors(slong* nb_acceptable_factors,
				      fmpz_poly_struct* acceptable_factors,
				      const fmpz_poly_t poly, slong ell)
{
  slong nb_factors = 0;
  slong max_nb_factors = fmpz_poly_degree(poly);
  fmpz_poly_struct* factors;
  slong* exps;
  fmpz_poly_t fac;
  slong d;
  slong k;

  factors = flint_malloc(max_nb_factors * sizeof(fmpz_poly_struct));
  for (k = 0; k < max_nb_factors; k++) fmpz_poly_init(&factors[k]);
  exps = flint_malloc(max_nb_factors * sizeof(slong));
  fmpz_poly_init(fac);
  
  modeq_factor_Q(&nb_factors, factors, exps, poly);
  *nb_acceptable_factors = 0;
  
  /* Isolate all factors of degree not 1 and dividing l+1 */
  for (k = 0; k < nb_factors; k++)
    {
      d = fmpz_poly_degree(&factors[k]);
      if (d > 1 && (ell+1) % d == 0) /* Degree of intermediate roots must divide l+1 */
	{
	  fmpz_poly_set(fac, &factors[k]);
	  fmpz_poly_primitive_part(fac, fac);
	  fmpz_poly_set(&acceptable_factors[*nb_acceptable_factors], fac);
	  nb_acceptable_factors++;
	  /* flint_printf("Degree %wd factor:\n", d);
	     fmpz_poly_print_pretty(fac, "x");
	     flint_printf("\n"); */
	}
    }

  for (k = 0; k < max_nb_factors; k++) fmpz_poly_clear(&factors[k]);
  flint_free(factors);
  flint_free(exps);
  fmpz_poly_clear(fac);
}
				      

/* Given a primitive, irreducible factor fac of &num_vec[0], compute
   an integer rescale with the following property: for any root j1 of
   fac, completed to (j1,j2,j3) following num_vec, the result of
   siegel_modeq_eval_C multiplied by rescale has integer coefficients
   in the number field defined by fac. */
static void
siegel_modeq_2step_rescale(fmpz_t rescale, const fmpz_poly_t fac,
			   const fmpz_poly_struct* num_vec, slong ell)
{
  fmpz* denoms;
  fmpq_poly_struct j_nf[3];
  fmpq_t coeff;
  fmpq_poly_t Q, R, fac_fmpq;
  fmpq_mat_t mult;
  fmpq_poly_t charpoly;
  mpz_t num1, den1;
  fmpz_t den;

  slong d = fmpz_poly_degree(fac);
  slong exp = (10*siegel_nb_cosets(ell))/3;
  slong k, m, n;

  denoms = _fmpz_vec_init(3);
  for (k = 0; k < 3; k++) fmpq_poly_init(&j_nf[3]);
  fmpq_init(coeff);
  fmpq_poly_init(Q);
  fmpq_poly_init(R);
  fmpq_poly_init(fac_fmpq);
  fmpq_mat_init(mult, d, d);
  fmpq_poly_init(charpoly);
  mpz_init(num1);
  mpz_init(den1);
  fmpz_init(den);

  modeq_isog_invariants_nf(j_nf, num_vec, 3, fac);
  fmpq_poly_set_fmpz_poly(fac_fmpq, fac);
  
  for (k = 0; k < 3; k++) fmpz_one(&denoms[k]);
  for (k = 0; k < 3; k++)
    {
      /* Get denominator of each j by computing charpoly of multiplication matrix */
      for (m = 0; m < d; m++)
	{
	  /* Copy coefficients of j_nf[k]*/
	  for (n = 0; n < d; n++)
	    {
	      fmpq_poly_get_coeff_fmpq(coeff, &j_nf[k], n);
	      fmpq_set(fmpq_mat_entry(mult, n, m), coeff);
	    }
	  /* Multiply by x and reduce mod fac */
	  fmpq_poly_shift_left(&j_nf[k], &j_nf[k], 1);
	  fmpq_poly_divrem(Q, R, &j_nf[k], fac_fmpq);
	}
      fmpq_mat_charpoly(charpoly, mult);
      for (m = 1; m <= d; m++)
	{
	  /* Adjust denominator and charpoly so that all coefficients
	     from x^d to x^(d-m) are integers */
	  fmpq_poly_get_coeff_fmpq(coeff, charpoly, d-m);
	  fmpq_get_mpz_frac(num1, den1, coeff); /* Should be in irreducible form */
	  fmpz_set_mpz(den, den1);
	  fmpz_mul(&denoms[k], &denoms[k], den);
	  /* Rescale charpoly accordingly */
	  fmpq_one(coeff);
	  fmpq_mul_fmpz(coeff, coeff, den);
	  fmpq_poly_rescale(charpoly, charpoly, coeff);
	}
      /* Check denominator is now 1 */
      fmpq_poly_get_denominator(den, charpoly);
      assert (fmpz_is_one(den));
    }

  /* Denominators are now computed; compute suitable powers */
  fmpz_lcm(&denoms[0], &denoms[0], &denoms[1]);
  fmpz_pow_ui(&denoms[0], &denoms[0], exp);
  fmpz_pow_ui(&denoms[2], &denoms[2], exp/2); /* Use the fact that degree in j3 is twice less */
  fmpz_lcm(rescale, &denoms[0], &denoms[2]);

  _fmpz_vec_clear(denoms, 3);
  for (k = 0; k < 3; k++) fmpq_poly_clear(&j_nf[3]);
  fmpq_clear(coeff);
  fmpq_poly_clear(Q);
  fmpq_poly_clear(R);
  fmpq_poly_clear(fac_fmpq);
  fmpq_mat_clear(mult);
  fmpq_poly_clear(charpoly);
  mpz_clear(num1);
  mpz_clear(den1);
  fmpz_clear(den);
}

/* Given a concatenation of d tuples of three complex polynomials
(Siegel modular equations of level ell), rescale them all, compute
trace of each of the three entries, and attempt to recognize a
polynomial with integer coefficients */
static int
siegel_modeq_2step_trace(fmpz_poly_struct* trace_nums, acb_poly_struct* nums_acb,
			 const fmpz_t rescale, slong d, slong ell, slong prec)
{
  slong n, i;
  acb_poly_t trace_num;
  acb_t rescale_acb;
  int success = 1;
  arf_t max_radius;

  acb_poly_init(trace_num);
  arf_init(max_radius);
  acb_init(rescale_acb);

  acb_set_fmpz(rescale_acb, rescale);  
  for (n = 0; n < 3 && success; n++)
    {
      acb_poly_zero(trace_num);
      for (i = 0; i < d; i++)
	{
	  acb_poly_scalar_mul(&nums_acb[n+3*i], &nums_acb[n+3*i], rescale_acb, prec);
	  acb_poly_add(trace_num, trace_num, &nums_acb[n+3*i], prec);
	}
      success = modeq_round_poly(&trace_nums[n], max_radius, trace_num,
				 siegel_nb_cosets(ell));
      if (!success)
	{ 
	  flint_printf("(siegel_modeq_2step_isog_invariants) Could not recognize integer coefficients: max radius ");
	  arf_printd(max_radius, 10);
	  /*flint_printf("\nConstant coefficient:\n");
	    arb_printd(acb_realref(acb_poly_get_coeff_ptr(trace_num, 0)), 100);*/
	  flint_printf("\n");
	}
    }
  return success;

  acb_poly_clear(trace_num);
  arf_clear(max_radius);
  acb_clear(rescale_acb);
}

/* Given a complex polynomial poly whose coefficients are integers in
   the number field defined by fac (for one of its complex
   embeddings), and a putative rational root of poly, check whether it
   is indeed a root */
static int
siegel_modeq_2step_confirm_root(const fmpq_t root, const fmpz_poly_t fac,
				const acb_poly_t poly, slong prec)
{
  return 1;
}

/* Attempt to compute rational 2-step l-isogenous invariants
corresponding to factor fac of the Siegel modular equations at current
precision. Return 1 if we should stop increasing precision, and set
success according to how the computation went. */
static int
siegel_modeq_2step_attempt(int* success, slong* nb_isog_j, fmpq* isog_j, const fmpq_t j,
			   const fmpz_poly_t fac, const fmpz_poly_struct* num_vec,
			   slong ell, slong prec)
{
  slong d;
  acb_poly_t fac_acb;
  slong max_nb_roots = siegel_nb_cosets(ell);
  slong nb_cplx_roots;
  acb_ptr roots;
  acb_ptr igusa_tuples;
  acb_poly_struct* nums_acb;
  acb_ptr dens_acb;
  fmpz_t rescale;
  fmpz_poly_struct trace_nums[3];
  slong nb_2step_roots;
  slong* exps;
  fmpq* ratl_roots;
  slong k;
  int stop;

  acb_poly_init(fac_acb);
  roots = flint_malloc(max_nb_roots * sizeof(acb_struct));
  for (k = 0; k < max_nb_roots; k++) acb_init(&roots[k]);
  igusa_tuples = flint_malloc(3 * max_nb_roots * sizeof(acb_struct));
  for (k = 0; k < 3*max_nb_roots; k++) acb_init(&igusa_tuples[k]);
  nums_acb = flint_malloc(3 * max_nb_roots * sizeof(acb_poly_struct));
  for (k = 0; k < 3*max_nb_roots; k++) acb_poly_init(&nums_acb[k]);
  dens_acb = flint_malloc(max_nb_roots * sizeof(acb_struct));
  for (k = 0; k < max_nb_roots; k++) acb_init(&dens_acb[k]);
  fmpz_init(rescale);
  for (k = 0; k < 3; k++) fmpz_poly_init(&trace_nums[k]);
  exps = flint_malloc(max_nb_roots * sizeof(slong));
  ratl_roots = flint_malloc(max_nb_roots * sizeof(fmpq));
  for (k = 0; k < max_nb_roots; k++) fmpq_init(&ratl_roots[k]);
  
  *success = 1;
  acb_poly_set_fmpz_poly(fac_acb, fac, prec);
  d = fmpz_poly_degree(fac);
  
  /* Try computing complex roots at current precision */
  nb_cplx_roots = acb_poly_find_roots(roots, fac_acb, NULL, 0, prec);
  if (nb_cplx_roots < d)
    {
      flint_printf("(siegel_modeq_2step_isog_invariants_Q) Not enough precision to isolate complex roots\n");
      *success = 0;
      stop = 0;
    }

  /* Compute list of l-isogenous invariants over C; evaluate
     modular equations over C at all these complex points */
  for (k = 0; k < d; k++)
    {
      if (*success) *success = modeq_isog_invariants_C(&igusa_tuples[3*k], num_vec,
						     &roots[k], 3, prec);
      if (*success) *success = siegel_modeq_eval_C(&nums_acb[3*k], &dens_acb[k],
						 &igusa_tuples[3*k], ell, prec);
    }

  /* Compute trace over Q of rescaled modular equations */
  if (*success)
    {
      siegel_modeq_2step_rescale(rescale, fac, num_vec, ell);
      *success = siegel_modeq_2step_trace(trace_nums, nums_acb, rescale, d, ell, prec);
      if (fmpz_poly_is_zero(&trace_nums[0]))
	{
	  flint_printf("(siegel_modeq_2step_isog_invariants) Trace over Q of modular equations at second step is zero: this case is not implemented\n");
	  *success = 0;
	  stop = 1;
	}
    }
  
  /* We have succeeded in rounding trace of nums_acb into integer polynomials.
     Now compute roots */
  if (*success)
    {
      *nb_isog_j = 0;
      stop = 1; /* In any case, stop increasing complex precision */
      modeq_roots_Q(&nb_2step_roots, ratl_roots, exps, &trace_nums[0]);
      flint_printf("(siegel_modeq_2step_isog_invariants) Found %wd rational roots\n", nb_2step_roots);
      for (k = 0; k < nb_2step_roots; k++)
	{
	  /* Compute isogenous Igusa invariants and add them to current list */
	  if (!fmpq_equal(&ratl_roots[k], &j[0])
	      && siegel_modeq_2step_confirm_root(&ratl_roots[k], fac,
						 &nums_acb[0], prec))		      
	    {
	      *success = *success && modeq_isog_invariants_Q(&isog_j[3*(*nb_isog_j)],
							     trace_nums, &ratl_roots[k], 3);
	      (*nb_isog_j)++;
	    }					  
	}
    }
  
  acb_poly_clear(fac_acb);
  for (k = 0; k < max_nb_roots; k++) acb_clear(&roots[k]);
  flint_free(roots);
  for (k = 0; k < 3*max_nb_roots; k++) acb_clear(&igusa_tuples[k]);
  flint_free(igusa_tuples);
  for (k = 0; k < 3*max_nb_roots; k++) acb_poly_clear(&nums_acb[k]);
  flint_free(nums_acb);
  for (k = 0; k < max_nb_roots; k++) acb_clear(&dens_acb[k]);
  flint_free(dens_acb);
  fmpz_clear(rescale);
  for (k = 0; k < 3; k++) fmpz_poly_clear(&trace_nums[k]);
  flint_free(exps);
  for (k = 0; k < max_nb_roots; k++) fmpq_clear(&ratl_roots[k]);
  flint_free(ratl_roots);
  
  return stop;
}

int siegel_modeq_2step_isog_invariants_Q(slong* nb_roots, fmpq* all_isog_j,
					 fmpq* j, slong ell)
{
  fmpz_poly_struct num_vec[3];
  fmpz_t den;
  
  slong nb_factors = 0;
  slong max_nb_factors = siegel_nb_cosets(ell);
  fmpz_poly_struct* factors;
  fmpz_poly_t fac;
  slong prec;
  int success;
  slong nb_2step_roots = 0;
  fmpq* isog_j;
  slong k, i;
  int res, stop;
  
  /* Init everything */
  for (k = 0; k < 3; k++) fmpz_poly_init(&num_vec[k]);
  fmpz_init(den);
  factors = flint_malloc(max_nb_factors * sizeof(fmpz_poly_struct));
  for (k = 0; k < max_nb_factors; k++) fmpz_poly_init(&factors[k]);
  fmpz_poly_init(fac);
  isog_j = flint_malloc(3 * max_nb_factors * sizeof(fmpq));
  for (k = 0; k < 3*max_nb_factors; k++) fmpq_init(&isog_j[k]);
  
  
  /* Evaluate modular equations over QQ */
  siegel_modeq_eval_Q(num_vec, den, j, ell);
  siegel_modeq_2step_acceptable_factors(&nb_factors, factors, &num_vec[0], ell);
  flint_printf("(siegel_modeq_2step_isog_invariants_Q) Identified %wd factors of correct degrees\n", nb_factors);

  /* For each of them, try to see if there are 2-step isogenous rational invariants. */
  *nb_roots = 0;
  res = 1; /* If we hit a problem, set to zero and exit */
  for (k = 0; k < nb_factors && res; k++)
    {
      fmpz_poly_set(fac, &factors[k]);
      prec = siegel_modeq_startprec_fmpq(j, ell);
      stop = 0;
      while (!stop) /* Keep increasing prec */
	{
	  flint_printf("(siegel_modeq_2step_isog_invariants_Q) Factor number %wd: starting new run at precision %wd\n", k, prec);
	  stop = siegel_modeq_2step_attempt(&success, &nb_2step_roots, isog_j, j,
					    fac, num_vec, ell, prec);
	  if (stop && success)
	    {
	      /* Copy isog_j to all_isog_j */
	      for (i = 0; i < nb_2step_roots; i++)
		{
		  fmpq_set(&all_isog_j[*nb_roots+i], &isog_j[i]);
		}
	      *nb_roots += nb_2step_roots;
	    }
	  if (stop && !success)
	    {
	      /* Computation has failed */
	      res = 0;
	    }
	  if (!stop)
	    {
	      prec = siegel_modeq_nextprec(prec);
	      if (prec > MODEQ_MAX_PREC)
		{
		  flint_printf("(siegel_modeq_2step_isog_invariants_Q) Reached maximal allowed precision %wd, abandon.\n", MODEQ_MAX_PREC);
		  stop = 1;
		  res = 0;
		} 
	    }
	} /* End loop to increase precision for given factor */
    } /* End loop over factors */

  /* Clear everything */  
  for (k = 0; k < 3; k++) fmpz_poly_clear(&num_vec[k]);
  fmpz_clear(den);
  for (k = 0; k < max_nb_factors; k++) fmpz_poly_clear(&factors[k]);
  flint_free(factors);
  fmpz_poly_clear(fac);
  for (k = 0; k < 3*max_nb_factors; k++) fmpq_clear(&isog_j[k]);
  flint_free(isog_j);  
  return res;
}
