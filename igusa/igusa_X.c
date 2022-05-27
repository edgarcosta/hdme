
#include "igusa.h"

static fmpq* Y12(fmpq* X) {return &X[4];}
static fmpq* X16(fmpq* X) {return &X[5];}
static fmpq* X18(fmpq* X) {return &X[6];}
static fmpq* X24(fmpq* X) {return &X[7];}
static fmpq* X28(fmpq* X) {return &X[8];}
static fmpq* X30(fmpq* X) {return &X[9];}
static fmpq* X36(fmpq* X) {return &X[10];}
static fmpq* X40(fmpq* X) {return &X[11];}
static fmpq* X42(fmpq* X) {return &X[12];}
static fmpq* X48(fmpq* X) {return &X[13];}

static void fmpq_div_si(fmpq_t r, fmpq_t x, slong c)
{
  fmpq_t t;
  fmpq_init(t);
  fmpq_set_si(t, c, 1);
  fmpq_div(r, x, t);
  fmpq_clear(t);  
}

static void igusa_X_q(fmpq* X, fmpq* I)
{
  fmpq_t c;
  fmpq_t temp;
  slong k;

  fmpq_init(c);
  fmpq_init(temp);

  for (k = 0; k < 4; k++) fmpq_set(&X[k], &I[k]);
  
  /* Y12 */
  fmpq_pow_si(c, igusa_psi4(I), 3);
  fmpq_pow_si(temp, igusa_psi6(I), 2);
  fmpq_sub(c, c, temp);
  fmpq_div_si(c, c, n_pow(2,6)*n_pow(3,3));
  fmpq_mul_si(temp, igusa_chi12(I), n_pow(2,4)*n_pow(3,2));
  fmpq_add(c, c, temp);
  fmpq_set(Y12(X), c);

  /* X16 */
  fmpq_mul(c, igusa_psi4(I), igusa_chi12(I));
  fmpq_mul(temp, igusa_psi6(I), igusa_chi10(I));
  fmpq_sub(c, c, temp);
  fmpq_div_si(c, c, 12);
  fmpq_set(X16(X), c);

  /* X18 */
  fmpq_mul(c, igusa_psi6(I), igusa_chi12(I));
  fmpq_mul(temp, igusa_psi4(I), igusa_psi4(I));
  fmpq_mul(temp, temp, igusa_chi10(I));
  fmpq_sub(c, c, temp);
  fmpq_div_si(c, c, 12);
  fmpq_set(X18(X), c);
  
  /* X24 */
  fmpq_mul(c, igusa_chi12(I), igusa_chi12(I));
  fmpq_mul(temp, igusa_chi10(I), igusa_chi10(I));
  fmpq_mul(temp, temp, igusa_psi4(I));
  fmpq_sub(c, c, temp);
  fmpq_div_si(c, c, 24);
  fmpq_set(X24(X), c);
  
  /* X28 */
  fmpq_mul(c, X24(X), igusa_psi4(I));
  fmpq_mul(temp, igusa_chi10(I), X18(X));
  fmpq_sub(c, c, temp);
  fmpq_div_si(c, c, 6);
  fmpq_set(X28(X), c);
  
  /* X30 */
  fmpq_mul(c, igusa_psi6(I), X24(X));
  fmpq_mul(temp, igusa_psi4(I), igusa_chi10(I));
  fmpq_mul(temp, temp, X16(X));
  fmpq_sub(c, c, temp);
  fmpq_div_si(c, c, 6);
  fmpq_set(X30(X), c);
  
  /* X36 */
  fmpq_mul(c, igusa_chi12(I), X24(X));
  fmpq_mul(temp, igusa_chi10(I), igusa_chi10(I));
  fmpq_mul(temp, temp, X16(X));
  fmpq_sub(c, c, temp);
  fmpq_div_si(c, c, 18);
  fmpq_div_si(X36(X), X36(X), 12);
  
  /* X40 */
  fmpq_mul(c, igusa_psi4(I), X36(X));
  fmpq_mul(temp, igusa_chi10(I), X30(X));
  fmpq_sub(c, c, temp);
  fmpq_div_si(c, c, 4);
  fmpq_set(X40(X), c);
  
  /* X42 */
  fmpq_mul(c, igusa_chi12(I), X30(X));
  fmpq_mul(temp, igusa_psi4(I), igusa_chi10(I));
  fmpq_mul(temp, temp, X28(X));
  fmpq_sub(c, c, temp);
  fmpq_div_si(c, c, 12);
  fmpq_set(X42(X), c);
  
  /* X48 */
  fmpq_mul(c, igusa_chi12(I), X36(X));
  fmpq_mul(temp, X24(X), X24(X));
  fmpq_sub(c, c, temp);
  fmpq_div_si(c, c, 4);
  fmpq_set(X48(X), c);

  fmpq_clear(c);
  fmpq_clear(temp);
}

void igusa_X(fmpq* X, fmpz* I)
{
  fmpq* J;
  fmpz_t one;
  slong k;

  fmpz_init(one);
  J = _fmpq_vec_init(4);
  
  fmpz_one(one);
  for (k = 0; k < 4; k++) fmpq_set_fmpz_frac(&J[k], &I[k], one);
  igusa_X_q(X, J);

  fmpz_clear(one);
  _fmpq_vec_clear(J, 4);
}
