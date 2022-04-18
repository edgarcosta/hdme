
#include "igusa.h"

void cov_from_igusa_fmpz(fmpz* I, fmpq* j)
{
  /* j1 = I4 I6'/I10; j2 = I2 I4^2/I10; j3 = I4^5/I10^2 I6 will be
   fixed by j1. I2 will be fixed by j2.  Solve j3 = I4^5/I10^2 by
   I4=n3*d3, I10=n3^2*d3^3 Assume j3 is nonzero. */
  mpz_t n, d;
  fmpz_t num, den, scal;
  
  mpz_init(n);
  mpz_init(d);
  fmpz_init(num);
  fmpz_init(den);
  fmpz_init(scal);
  
  fmpq_get_mpz_frac(n, d, &j[2]);
  fmpz_set_mpz(num, n);
  fmpz_set_mpz(den, d);
  fmpz_mul(&I[1], num, den);
  fmpz_mul(&I[3], &I[1], &I[1]);
  fmpz_mul(&I[3], &I[3], den);
  
  /* Prepare others: I6' <- n3*d3^2, I2 <- d3 */
  fmpz_set(&I[0], den);
  fmpz_mul(&I[2], num, den);
  fmpz_mul(&I[2], &I[2], den);

  /* Adjust I2 */
  fmpq_get_mpz_frac(n, d, &j[1]);
  fmpz_set_mpz(num, n);
  fmpz_set_mpz(den, d);
  fmpz_mul(&I[0], &I[0], num);
  /* We want to divide I2 by den; if not possible, rescale */
  if (!fmpz_divisible(&I[0], den)) cov_rescale_fmpz(I, I, den);
  fmpz_divexact(&I[0], &I[0], den);

  /* Adjust I6' */
  fmpq_get_mpz_frac(n, d, &j[0]);
  fmpz_set_mpz(num, n);
  fmpz_set_mpz(den, d);
  fmpz_mul(&I[2], &I[2], num);
  /* We want to divide I6' by den */
  if (!fmpz_divisible(&I[2], den)) cov_rescale_fmpz(I, I, den);
  fmpz_divexact(&I[2], &I[2], den);

  cov_normalize_fmpz(I, I);
  
  mpz_clear(n);
  mpz_clear(d);
  fmpz_clear(num);
  fmpz_clear(den);
  fmpz_clear(scal);
}
