
#include "modular.h"

static void set_pol(fmpz_mpoly_t pol, const modeq_ctx_t ctx, slong j)
{
  slong k = 0;
  fmpz_mpoly_t term;

  fmpz_mpoly_init(term, modeq_ctx_ctx(ctx));
  fmpz_mpoly_set(pol, modeq_ctx_monomial(ctx, 0), modeq_ctx_ctx(ctx));
  
  while (j > 0)
    {
      fmpz_mpoly_set(term, modeq_ctx_monomial(ctx, k+1), modeq_ctx_ctx(ctx));
      if (j%2 == 1) fmpz_mpoly_add(pol, pol, term, modeq_ctx_ctx(ctx));
      j = j/2;
      k++;
    }
  
  fmpz_mpoly_clear(term, modeq_ctx_ctx(ctx));
}

int modeq_ctx_choose(modeq_ctx_t ctx, acb_srcptr I, slong nb, slong prec)
{
  slong j, k, l, e;
  acb_srcptr ptr; /* Do not initialize */
  slong weights[4] = IGUSA_HALFWEIGHTS;
  slong exps[4];
  fmpz_mpoly_t num, den;
  acb_ptr evden;
  acb_ptr evnum;
  int res;
  slong wt = 60;
  int v = MODEQ_VERBOSE;

  fmpz_mpoly_init(num, modeq_ctx_ctx(ctx));
  fmpz_mpoly_init(den, modeq_ctx_ctx(ctx));
  evden = _acb_vec_init(nb);
  evnum = _acb_vec_init(nb);

  /* Do we have: I4 is never zero? */
  res = 1;
  for (j = 0; j < nb; j++)
    {
      if (acb_contains_zero(igusa_I4(&I[4*j])))
	{
	  res = 0;
	  break;
	}
    }
  if (res) wt = 20;

  /* Do we have: I4=I6=0 and I6prime=I10=0 never happen? */
  res = 1;
  for (j = 0; j < nb; j++)
    {
      ptr = &I[4*j];
      if (acb_contains_zero(igusa_I4(ptr))
	  && acb_contains_zero(igusa_I6prime(ptr)))
	{
	  res = 0;
	  break;
	}
      if (acb_contains_zero(igusa_I6prime(ptr))
	  && acb_contains_zero(igusa_I10(ptr)))
	{
	  res = 0;
	  break;
	}
    }
  if (res && wt == 60) wt = 30;

  /* Set monomials in ctx */
  modeq_ctx_weight(ctx) = wt;
  modeq_ctx_nb(ctx) = igusa_nb_base_monomials(wt);
  for (k = 0; k < modeq_ctx_nb(ctx); k++)
    {
      igusa_base_exps(exps, wt, k);
      igusa_base_monomial(modeq_ctx_monomial(ctx, k), wt, k,
			  modeq_ctx_ctx(ctx));
    }

  /* Collect possibly equal pairs */
  for (j = 0; j < nb; j++)
    {
      for (k = j+1; k < nb; k++)
	{
	  if (!cov_distinct(&I[4*j], &I[4*k], 4, weights, prec))
	    {
	      if (v)
		{
		  flint_printf("(modeq_ctx_choose) Found suspicious pair: %wd, %wd\n",
			       j, k);
		}
	      modeq_ctx_add_pair(ctx, j, k);
	    }
	}
    }

  /* Choose denominator: should never vanish */
  /* Each vanishing, for ell fixed, defines a hypersurface in
     the moduli space. Choosing 8 of them should be sufficient. */
  e = 3;
  for (j = 0; j < n_pow(2,e); j++)
    {
      set_pol(den, ctx, j);
      res = 1;
      for (k = 0; k < nb; k++)
	{
	  cov_mpoly_eval(&evden[k], den, &I[4*k], modeq_ctx_ctx(ctx), prec);
	  if (acb_contains_zero(&evden[k]))
	    {
	      res = 0;
	      break;
	    }
	}
      if (res)
	{
	  fmpz_mpoly_set(modeq_ctx_den(ctx), den, modeq_ctx_ctx(ctx));
	  break;
	}
    }
  if (!res && v)
    {
      flint_printf("(modeq_ctx_choose) Warning: could not find suitable denominator\n");
    }
    
  /* Choose numerator: num/den should be distinct on non-"equal" pairs */
  /* Again, 8 of them should be sufficient. */
  if (res)
    {
      for (j = 0; j < n_pow(2, e); j++)
	{
	  set_pol(num, ctx, j);
	  for (k = 0; k < nb; k++)
	    {
	      cov_mpoly_eval(&evnum[k], num, &I[4*k], modeq_ctx_ctx(ctx), prec);
	      acb_div(&evnum[k], &evnum[k], &evden[k], prec);
	    }
	  res = 1;
	  for (k = 0; k < nb; k++)
	    {
	      for (l = k+1; l < nb; l++)
		{
		  /* For performance, use that pairs are ordered? */
		  if (acb_overlaps(&evnum[k], &evnum[l])
		      && !modeq_ctx_is_pair(k, l, ctx))
		    {
		      res = 0;
		      break;
		    }
		}
	      if (!res) break;
	    }
	  if (res) /* Numerator no. j works */
	    {
	      fmpz_mpoly_set(modeq_ctx_num(ctx), num, modeq_ctx_ctx(ctx));
	      break;
	    }
	}
      if (!res && v)
	{
	  flint_printf("(modeq_ctx_choose) Warning: could not find suitable numerator\n");
	}      
    }
    
  fmpz_mpoly_clear(num, modeq_ctx_ctx(ctx));
  fmpz_mpoly_clear(den, modeq_ctx_ctx(ctx));
  _acb_vec_clear(evden, nb);
  _acb_vec_clear(evnum, nb);
  return res;
}
