
#ifndef SIEGEL_H
#define SIEGEL_H

#include <stdio.h>
#include "flint/fmpz_mat.h"
#include "flint/ulong_extras.h"
#include "acb_mat.h"
#include "arb_mat.h"
#include "acb_mat_extras.h"

typedef struct
{
  fmpz_mat_struct a;
  fmpz_mat_struct b;
  fmpz_mat_struct c;
  fmpz_mat_struct d;
  slong g;
} sp2gz_struct;

typedef sp2gz_struct sp2gz_t[1];

ARB_INLINE void
sp2gz_init(sp2gz_t m, slong g)
{
  fmpz_mat_init(&m->a, g, g);
  fmpz_mat_init(&m->b, g, g);
  fmpz_mat_init(&m->c, g, g);
  fmpz_mat_init(&m->d, g, g);
  fmpz_mat_one(&m->a);
  fmpz_mat_one(&m->d);
  m->g = g;
}

ARB_INLINE void
sp2gz_clear(sp2gz_t m)
{
  fmpz_mat_clear(&m->a);
  fmpz_mat_clear(&m->b);
  fmpz_mat_clear(&m->c);
  fmpz_mat_clear(&m->d);
}

ARB_INLINE void
sp2gz_swap(sp2gz_t m, sp2gz_t n)
{
  sp2gz_struct p = *m;
  *m = *n;
  *n = p;
}

ARB_INLINE void
sp2gz_set(sp2gz_t m, const sp2gz_t n)
{
  fmpz_mat_set(&m->a, &n->a);
  fmpz_mat_set(&m->b, &n->b);
  fmpz_mat_set(&m->c, &n->c);
  fmpz_mat_set(&m->d, &n->d);
}

void sp2gz_set_mat(sp2gz_t m, const fmpz_mat_t n);

void sp2gz_get_mat(fmpz_mat_t m, const sp2gz_t n);

ARB_INLINE void
sp2gz_one(sp2gz_t m)
{
  fmpz_mat_one(&m->a);
  fmpz_mat_zero(&m->b);
  fmpz_mat_zero(&m->c);
  fmpz_mat_one(&m->d);
}

ARB_INLINE void
sp2gz_J(sp2gz_t m)
{
  fmpz_mat_zero(&m->a);
  fmpz_mat_one(&m->b);
  fmpz_mat_one(&m->c);
  fmpz_mat_neg(&m->c, &m->c);
  fmpz_mat_zero(&m->d);
}

void sp2gz_fprint(FILE * file, const sp2gz_t m);

ARB_INLINE void
sp2gz_print(const sp2gz_t m)
{
  sp2gz_fprint(stdout, m);
}

ARB_INLINE int
sp2gz_equal(const sp2gz_t m, const sp2gz_t n)
{
  return fmpz_mat_equal(&m->a, &n->a)
    &&   fmpz_mat_equal(&m->b, &n->b)
    &&   fmpz_mat_equal(&m->c, &n->c)
    &&   fmpz_mat_equal(&m->d, &n->d);
}

void sp2gz_mul(sp2gz_t r, const sp2gz_t m, const sp2gz_t n);

void sp2gz_inv(sp2gz_t r, const sp2gz_t m);

int sp2gz_is_one(const sp2gz_t m);

int sp2gz_is_J(const sp2gz_t m);

int sp2gz_is_correct(const sp2gz_t m);

void sp2gz_set_diagonal(sp2gz_t m, const fmpz_mat_t u);

void sp2gz_randtest_triangular(sp2gz_t m, flint_rand_t state, slong bits);

void sp2gz_randtest_diagonal(sp2gz_t m, flint_rand_t state, slong bits);

void sp2gz_randtest(sp2gz_t m, flint_rand_t state, slong bits);

void siegel_halfspace_randtest(acb_mat_t z, flint_rand_t state, slong prec);

void siegel_star(acb_mat_t w, const sp2gz_t m, const acb_mat_t z, slong prec);

int siegel_transform(acb_mat_t w, const sp2gz_t m, const acb_mat_t z, slong prec);

int siegel_is_real_reduced(const acb_mat_t z, const arb_t tol, slong prec);

int siegel_not_real_reduced(const acb_mat_t z, slong prec);

int siegel_reduce_real(acb_mat_t w, sp2gz_t u, const acb_mat_t z,
		       const arb_t tol, slong prec);

slong siegel_nb_test_matrices(slong g);

void siegel_test_matrix(sp2gz_t u, slong j);

int siegel_fundamental_domain(acb_mat_t w, sp2gz_t m,
			      const acb_mat_t z, const arb_t tol, slong prec);

int siegel_is_in_fundamental_domain(const acb_mat_t z, const arb_t tol, slong prec);

int siegel_not_in_fundamental_domain(const acb_mat_t z, slong prec);

int siegel_is_weakly_reduced(const acb_mat_t z, const arb_t tol, slong prec);

void siegel_fundamental_domain_randtest(acb_mat_t z, flint_rand_t state, slong prec);


#endif 
