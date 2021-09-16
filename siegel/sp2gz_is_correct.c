
#include "siegel.h"

int
sp2gz_is_correct(const sp2gz_t m)
{
  fmpz_mat_t prod1, prod2;
  int eq1, eq2, eq3;

  fmpz_mat_init(prod1, m->g, m->g);
  fmpz_mat_init(prod2, m->g, m->g);

  fmpz_mat_transpose(prod1, &m->a);
  fmpz_mat_mul(prod1, prod1, &m->c);
  fmpz_mat_transpose(prod2, &m->c);
  fmpz_mat_mul(prod2, prod2, &m->a);
  fmpz_mat_sub(prod1, prod1, prod2);
  eq1 = fmpz_mat_is_zero(prod1);

  fmpz_mat_transpose(prod1, &m->b);
  fmpz_mat_mul(prod1, prod1, &m->d);
  fmpz_mat_transpose(prod2, &m->d);
  fmpz_mat_mul(prod2, prod2, &m->b);
  fmpz_mat_sub(prod1, prod1, prod2);
  eq2 = fmpz_mat_is_zero(prod1);

  fmpz_mat_transpose(prod1, &m->a);
  fmpz_mat_mul(prod1, prod1, &m->d);
  fmpz_mat_transpose(prod2, &m->c);
  fmpz_mat_mul(prod2, prod2, &m->b);
  fmpz_mat_sub(prod1, prod1, prod2);
  eq3 = fmpz_mat_is_one(prod1);

  fmpz_mat_clear(prod1);
  fmpz_mat_clear(prod2);
  
  return eq1 && eq2 && eq3;
}
