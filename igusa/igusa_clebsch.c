
#include "igusa.h"

void igusa_clebsch(acb_ptr ABCD, acb_srcptr I, slong prec)
{
  acb_t I2, I4, I6, I10, A, B, C, D, temp;
  acb_ptr A_pows;
  acb_ptr B_pows;
  acb_ptr C_pows;
  acb_ptr I2_pows;
  acb_ptr I4_pows;
  acb_ptr I6_pows;
  acb_ptr I10_pows;

  acb_init(I2);
  acb_init(I4);
  acb_init(I6);
  acb_init(I10);
  acb_init(A);
  acb_init(B);
  acb_init(C);
  acb_init(D);
  acb_init(temp);
  
  acb_set(I2, &I[0]);
  acb_set(I4, &I[1]);
  acb_set(I10, &I[3]);
  igusa_I6(I6, I, prec);

  /* Auto-generated code */

  
  /* Init all power lists */
  I2_pows = _acb_vec_init(2);

  /* Precompute all power lists */
  acb_one(&I2_pows[0]);
  acb_mul(&I2_pows[1], &I2_pows[0], I2, prec);

  /* Add successive terms to res */
  acb_zero(A);

  acb_one(temp);
  acb_mul(temp, temp, &I2_pows[1], prec);
  acb_mul_si(temp, temp, -1, prec);
  acb_div_si(temp, temp, 120, prec);
  acb_add(A, A, temp, prec);

  /* Clear all power lists */
  _acb_vec_clear(I2_pows, 2);


  /* Init all power lists */
  A_pows = _acb_vec_init(3);
  I4_pows = _acb_vec_init(2);

  /* Precompute all power lists */
  acb_one(&A_pows[0]);
  acb_mul(&A_pows[1], &A_pows[0], A, prec);
  acb_mul(&A_pows[2], &A_pows[1], A, prec);

  acb_one(&I4_pows[0]);
  acb_mul(&I4_pows[1], &I4_pows[0], I4, prec);

  /* Add successive terms to res */
  acb_zero(B);

  acb_one(temp);
  acb_mul(temp, temp, &A_pows[2], prec);
  acb_mul_si(temp, temp, 8, prec);
  acb_div_si(temp, temp, 75, prec);
  acb_add(B, B, temp, prec);

  acb_one(temp);
  acb_mul(temp, temp, &I4_pows[1], prec);
  acb_mul_si(temp, temp, 1, prec);
  acb_div_si(temp, temp, 6750, prec);
  acb_add(B, B, temp, prec);

  /* Clear all power lists */
  _acb_vec_clear(A_pows, 3);
  _acb_vec_clear(I4_pows, 2);


  /* Init all power lists */
  A_pows = _acb_vec_init(4);
  B_pows = _acb_vec_init(2);
  I6_pows = _acb_vec_init(2);

  /* Precompute all power lists */
  acb_one(&A_pows[0]);
  acb_mul(&A_pows[1], &A_pows[0], A, prec);
  acb_mul(&A_pows[2], &A_pows[1], A, prec);
  acb_mul(&A_pows[3], &A_pows[2], A, prec);

  acb_one(&B_pows[0]);
  acb_mul(&B_pows[1], &B_pows[0], B, prec);

  acb_one(&I6_pows[0]);
  acb_mul(&I6_pows[1], &I6_pows[0], I6, prec);

  /* Add successive terms to res */
  acb_zero(C);

  acb_one(temp);
  acb_mul(temp, temp, &A_pows[3], prec);
  acb_mul_si(temp, temp, -16, prec);
  acb_div_si(temp, temp, 375, prec);
  acb_add(C, C, temp, prec);

  acb_one(temp);
  acb_mul(temp, temp, &A_pows[1], prec);
  acb_mul(temp, temp, &B_pows[1], prec);
  acb_mul_si(temp, temp, 8, prec);
  acb_div_si(temp, temp, 15, prec);
  acb_add(C, C, temp, prec);

  acb_one(temp);
  acb_mul(temp, temp, &I6_pows[1], prec);
  acb_mul_si(temp, temp, 1, prec);
  acb_div_si(temp, temp, 202500, prec);
  acb_add(C, C, temp, prec);

  /* Clear all power lists */
  _acb_vec_clear(A_pows, 4);
  _acb_vec_clear(B_pows, 2);
  _acb_vec_clear(I6_pows, 2);


  /* Init all power lists */
  A_pows = _acb_vec_init(6);
  B_pows = _acb_vec_init(3);
  C_pows = _acb_vec_init(2);
  I10_pows = _acb_vec_init(2);

  /* Precompute all power lists */
  acb_one(&A_pows[0]);
  acb_mul(&A_pows[1], &A_pows[0], A, prec);
  acb_mul(&A_pows[2], &A_pows[1], A, prec);
  acb_mul(&A_pows[3], &A_pows[2], A, prec);
  acb_mul(&A_pows[4], &A_pows[3], A, prec);
  acb_mul(&A_pows[5], &A_pows[4], A, prec);

  acb_one(&B_pows[0]);
  acb_mul(&B_pows[1], &B_pows[0], B, prec);
  acb_mul(&B_pows[2], &B_pows[1], B, prec);

  acb_one(&C_pows[0]);
  acb_mul(&C_pows[1], &C_pows[0], C, prec);

  acb_one(&I10_pows[0]);
  acb_mul(&I10_pows[1], &I10_pows[0], I10, prec);

  /* Add successive terms to res */
  acb_zero(D);

  acb_one(temp);
  acb_mul(temp, temp, &A_pows[5], prec);
  acb_mul_si(temp, temp, -128, prec);
  acb_div_si(temp, temp, 9375, prec);
  acb_add(D, D, temp, prec);

  acb_one(temp);
  acb_mul(temp, temp, &A_pows[3], prec);
  acb_mul(temp, temp, &B_pows[1], prec);
  acb_mul_si(temp, temp, 16, prec);
  acb_div_si(temp, temp, 75, prec);
  acb_add(D, D, temp, prec);

  acb_one(temp);
  acb_mul(temp, temp, &A_pows[1], prec);
  acb_mul(temp, temp, &B_pows[2], prec);
  acb_mul_si(temp, temp, -2, prec);
  acb_div_si(temp, temp, 3, prec);
  acb_add(D, D, temp, prec);

  acb_one(temp);
  acb_mul(temp, temp, &A_pows[2], prec);
  acb_mul(temp, temp, &C_pows[1], prec);
  acb_mul_si(temp, temp, 16, prec);
  acb_div_si(temp, temp, 45, prec);
  acb_add(D, D, temp, prec);

  acb_one(temp);
  acb_mul(temp, temp, &B_pows[1], prec);
  acb_mul(temp, temp, &C_pows[1], prec);
  acb_mul_si(temp, temp, -4, prec);
  acb_div_si(temp, temp, 3, prec);
  acb_add(D, D, temp, prec);

  acb_one(temp);
  acb_mul(temp, temp, &I10_pows[1], prec);
  acb_mul_si(temp, temp, -1, prec);
  acb_div_si(temp, temp, 4556250, prec);
  acb_add(D, D, temp, prec);

  /* Clear all power lists */
  _acb_vec_clear(A_pows, 6);
  _acb_vec_clear(B_pows, 3);
  _acb_vec_clear(C_pows, 2);
  _acb_vec_clear(I10_pows, 2);

  acb_set(&ABCD[0], A);
  acb_set(&ABCD[1], B);
  acb_set(&ABCD[2], C);
  acb_set(&ABCD[3], D);
  
  acb_clear(I2);
  acb_clear(I4);
  acb_clear(I6);
  acb_clear(I10);
  acb_clear(A);
  acb_clear(B);
  acb_clear(C);
  acb_clear(D);
  acb_clear(temp);
}
