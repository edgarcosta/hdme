
#include "hilbert.h"

void gundlach_from_hilbert_param(fmpq* res, const fmpq* mn, slong delta)
{
  fmpq_t g;
  fmpq_t h;
  fmpq_t temp, temp2;

  if (delta != 5)
    {
      flint_printf("(gundlach_from_hilbert_param) Gundlach invariants only implemented for delta = 5\n");
      fflush(stdout);
      flint_abort();
    }

  fmpq_init(g);
  fmpq_init(h);
  fmpq_init(temp);
  fmpq_init(temp2);

  fmpq_pow_si(g, &mn[0], 2);
  fmpq_pow_si(temp, &mn[1], 2);
  fmpq_mul_si(temp, temp, -5);
  fmpq_add(g, g, temp);
  fmpq_add_si(g, g, -9);
  fmpq_set_si(temp, 30, 1);
  fmpq_div(g, g, temp);

  fmpq_mul_si(h, &mn[0], 3);
  fmpq_mul_si(temp, g, 10);
  fmpq_add_si(temp, temp, 3);
  fmpq_mul(h, h, temp);
  fmpq_mul_si(temp, g, 15);
  fmpq_add_si(temp, temp, 2);
  fmpq_mul(h, h, temp);
  fmpq_pow_si(temp, g, 2);
  fmpq_mul_si(temp, temp, 250);
  fmpq_mul_si(temp2, g, 75);
  fmpq_add(temp, temp, temp2);
  fmpq_add_si(temp, temp, 6);
  fmpq_mul_si(temp, temp, 9);
  fmpq_add(h, h, temp);
  fmpq_set_si(temp, 6250, 1);
  fmpq_div(h, h, temp);

  fmpq_mul_si(g, g, -6);
  fmpq_pow_si(&res[0], g, 5);
  fmpq_pow_si(temp, h, 2);
  fmpq_div(&res[0], &res[0], temp);
  fmpq_pow_si(&res[1], g, 2);
  fmpq_div(&res[1], &res[1], h);

  fmpq_clear(g);
  fmpq_clear(h);
  fmpq_clear(temp);
  fmpq_clear(temp2);
}
