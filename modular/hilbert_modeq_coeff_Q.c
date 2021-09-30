
#include "hilbert.h"

static int
cont_frac_step(fmpz_t r, arf_t next, const arf_t current, slong prec)
{
  int res = 0;    
  arf_get_fmpz(r, current, ARF_RND_FLOOR);
  arf_sub_fmpz(next, current, r, prec, ARF_RND_NEAR);
  if (arf_cmp_2exp_si(next, -50) < 0)
    {
      res = 1;
    }
  else
    {
      arf_ui_div(next, 1, next, prec, ARF_RND_NEAR);
    }
  return res;
}

static slong
cont_frac_max_nb_steps(const acb_t x, slong prec)
{
  return prec/2;
}

static void
cont_frac_get_fmpq(fmpq_t c, fmpz* r_vec, slong nb_steps)
{
  slong k;
  fmpq_zero(c);
  fmpq_add_fmpz(c, c, &r_vec[nb_steps-1]);
  for (k = nb_steps-2; k >= 0; k--)
    {
      fmpq_inv(c, c);
      fmpq_add_fmpz(c, c, &r_vec[k]);
    }
}

int hilbert_modeq_coeff_Q(fmpq_t c, fmpz_t den, const acb_t x,
			  const fmpz_t probable_den, slong prec)
{
  acb_t z;
  arf_t current;
  mpz_t n, d;
  int res = 0;
  int stop = 0;
  slong max_steps = cont_frac_max_nb_steps(x, prec);
  fmpz* r_vec;
  slong k;

  acb_init(z);
  arf_init(current);
  mpz_init(n);
  mpz_init(d);
  r_vec = _fmpz_vec_init(max_steps);
  
  acb_mul_fmpz(z, x, probable_den, prec);
  
  if (!arb_contains_zero(acb_imagref(z)))
    {
      res = 0;
    }
  if (res)
    {
      arf_set(current, arb_midref(acb_realref(z)));
      k = 0;
      while (!stop && (k < max_steps))
	{
	  stop = cont_frac_step(&r_vec[k], current, current, prec);
	  k++;
	}
      if (k == max_steps)
	{
	  res = 0;
	}
      else
	{
	  res = 1;
	  cont_frac_get_fmpq(c, r_vec, k);
	  fmpq_div_fmpz(c, c, probable_den);
	  fmpq_get_mpz_frac(n, d, c);
	  fmpz_set_mpz(den, d);
	}
    }

  acb_clear(z);
  arf_clear(current);
  mpz_clear(n);
  mpz_clear(d);
  _fmpz_vec_clear(r_vec, max_steps);
  return res;
}
