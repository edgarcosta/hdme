
#include "igusa.h"

void mestre_conic(acb_ptr conic, acb_srcptr ABCD, const acb_t U, const acb_t I10, slong prec)
{
  acb_ptr A_coefs;
  acb_ptr C_coefs;
  acb_t temp;

  acb_ptr A_pows;
  acb_ptr B_pows;
  acb_ptr C_pows;
  acb_ptr D_pows;

  A_coefs = _acb_vec_init(6);
  C_coefs = _acb_vec_init(6);
  acb_init(temp);
  
  /* Init all power lists */
  A_pows = _acb_vec_init(2);
  B_pows = _acb_vec_init(2);
  C_pows = _acb_vec_init(2);

  /* Precompute all power lists */
  acb_one(&A_pows[0]);
  acb_mul(&A_pows[1], &A_pows[0], &ABCD[0], prec);

  acb_one(&B_pows[0]);
  acb_mul(&B_pows[1], &B_pows[0], &ABCD[1], prec);

  acb_one(&C_pows[0]);
  acb_mul(&C_pows[1], &C_pows[0], &ABCD[2], prec);

  /* Add successive terms to res */
  acb_zero(&A_coefs[0]);

  acb_one(temp);
  acb_mul(temp, temp, &A_pows[1], prec);
  acb_mul(temp, temp, &B_pows[1], prec);
  acb_mul_si(temp, temp, 1, prec);
  acb_div_si(temp, temp, 3, prec);
  acb_add(&A_coefs[0], &A_coefs[0], temp, prec);

  acb_one(temp);
  acb_mul(temp, temp, &C_pows[1], prec);
  acb_mul_si(temp, temp, 2, prec);
  acb_div_si(temp, temp, 1, prec);
  acb_add(&A_coefs[0], &A_coefs[0], temp, prec);

  /* Clear all power lists */
  _acb_vec_clear(A_pows, 2);
  _acb_vec_clear(B_pows, 2);
  _acb_vec_clear(C_pows, 2);


  /* Init all power lists */
  D_pows = _acb_vec_init(2);

  /* Precompute all power lists */
  acb_one(&D_pows[0]);
  acb_mul(&D_pows[1], &D_pows[0], &ABCD[3], prec);

  /* Add successive terms to res */
  acb_zero(&A_coefs[1]);

  acb_one(temp);
  acb_mul(temp, temp, &D_pows[1], prec);
  acb_mul_si(temp, temp, 1, prec);
  acb_div_si(temp, temp, 1, prec);
  acb_add(&A_coefs[1], &A_coefs[1], temp, prec);

  /* Clear all power lists */
  _acb_vec_clear(D_pows, 2);


  /* Init all power lists */
  A_pows = _acb_vec_init(2);
  B_pows = _acb_vec_init(3);
  C_pows = _acb_vec_init(3);
  D_pows = _acb_vec_init(2);

  /* Precompute all power lists */
  acb_one(&A_pows[0]);
  acb_mul(&A_pows[1], &A_pows[0], &ABCD[0], prec);

  acb_one(&B_pows[0]);
  acb_mul(&B_pows[1], &B_pows[0], &ABCD[1], prec);
  acb_mul(&B_pows[2], &B_pows[1], &ABCD[1], prec);

  acb_one(&C_pows[0]);
  acb_mul(&C_pows[1], &C_pows[0], &ABCD[2], prec);
  acb_mul(&C_pows[2], &C_pows[1], &ABCD[2], prec);

  acb_one(&D_pows[0]);
  acb_mul(&D_pows[1], &D_pows[0], &ABCD[3], prec);

  /* Add successive terms to res */
  acb_zero(&A_coefs[2]);

  acb_one(temp);
  acb_mul(temp, temp, &B_pows[2], prec);
  acb_mul(temp, temp, &C_pows[1], prec);
  acb_mul_si(temp, temp, 2, prec);
  acb_div_si(temp, temp, 9, prec);
  acb_add(&A_coefs[2], &A_coefs[2], temp, prec);

  acb_one(temp);
  acb_mul(temp, temp, &A_pows[1], prec);
  acb_mul(temp, temp, &C_pows[2], prec);
  acb_mul_si(temp, temp, 2, prec);
  acb_div_si(temp, temp, 9, prec);
  acb_add(&A_coefs[2], &A_coefs[2], temp, prec);

  acb_one(temp);
  acb_mul(temp, temp, &B_pows[1], prec);
  acb_mul(temp, temp, &D_pows[1], prec);
  acb_mul_si(temp, temp, 1, prec);
  acb_div_si(temp, temp, 2, prec);
  acb_add(&A_coefs[2], &A_coefs[2], temp, prec);

  /* Clear all power lists */
  _acb_vec_clear(A_pows, 2);
  _acb_vec_clear(B_pows, 3);
  _acb_vec_clear(C_pows, 3);
  _acb_vec_clear(D_pows, 2);


  /* Init all power lists */
  A_pows = _acb_vec_init(2);
  B_pows = _acb_vec_init(4);
  C_pows = _acb_vec_init(3);

  /* Precompute all power lists */
  acb_one(&A_pows[0]);
  acb_mul(&A_pows[1], &A_pows[0], &ABCD[0], prec);

  acb_one(&B_pows[0]);
  acb_mul(&B_pows[1], &B_pows[0], &ABCD[1], prec);
  acb_mul(&B_pows[2], &B_pows[1], &ABCD[1], prec);
  acb_mul(&B_pows[3], &B_pows[2], &ABCD[1], prec);

  acb_one(&C_pows[0]);
  acb_mul(&C_pows[1], &C_pows[0], &ABCD[2], prec);
  acb_mul(&C_pows[2], &C_pows[1], &ABCD[2], prec);

  /* Add successive terms to res */
  acb_zero(&A_coefs[3]);

  acb_one(temp);
  acb_mul(temp, temp, &B_pows[3], prec);
  acb_mul_si(temp, temp, 1, prec);
  acb_div_si(temp, temp, 3, prec);
  acb_add(&A_coefs[3], &A_coefs[3], temp, prec);

  acb_one(temp);
  acb_mul(temp, temp, &A_pows[1], prec);
  acb_mul(temp, temp, &B_pows[1], prec);
  acb_mul(temp, temp, &C_pows[1], prec);
  acb_mul_si(temp, temp, 4, prec);
  acb_div_si(temp, temp, 9, prec);
  acb_add(&A_coefs[3], &A_coefs[3], temp, prec);

  acb_one(temp);
  acb_mul(temp, temp, &C_pows[2], prec);
  acb_mul_si(temp, temp, 2, prec);
  acb_div_si(temp, temp, 3, prec);
  acb_add(&A_coefs[3], &A_coefs[3], temp, prec);

  /* Clear all power lists */
  _acb_vec_clear(A_pows, 2);
  _acb_vec_clear(B_pows, 4);
  _acb_vec_clear(C_pows, 3);


  /* Init all power lists */
  D_pows = _acb_vec_init(2);

  /* Precompute all power lists */
  acb_one(&D_pows[0]);
  acb_mul(&D_pows[1], &D_pows[0], &ABCD[3], prec);

  /* Add successive terms to res */
  acb_zero(&A_coefs[4]);

  acb_one(temp);
  acb_mul(temp, temp, &D_pows[1], prec);
  acb_mul_si(temp, temp, 1, prec);
  acb_div_si(temp, temp, 1, prec);
  acb_add(&A_coefs[4], &A_coefs[4], temp, prec);

  /* Clear all power lists */
  _acb_vec_clear(D_pows, 2);


  /* Init all power lists */
  A_pows = _acb_vec_init(2);
  B_pows = _acb_vec_init(3);
  C_pows = _acb_vec_init(2);

  /* Precompute all power lists */
  acb_one(&A_pows[0]);
  acb_mul(&A_pows[1], &A_pows[0], &ABCD[0], prec);

  acb_one(&B_pows[0]);
  acb_mul(&B_pows[1], &B_pows[0], &ABCD[1], prec);
  acb_mul(&B_pows[2], &B_pows[1], &ABCD[1], prec);

  acb_one(&C_pows[0]);
  acb_mul(&C_pows[1], &C_pows[0], &ABCD[2], prec);

  /* Add successive terms to res */
  acb_zero(&A_coefs[5]);

  acb_one(temp);
  acb_mul(temp, temp, &B_pows[2], prec);
  acb_mul_si(temp, temp, 2, prec);
  acb_div_si(temp, temp, 3, prec);
  acb_add(&A_coefs[5], &A_coefs[5], temp, prec);

  acb_one(temp);
  acb_mul(temp, temp, &A_pows[1], prec);
  acb_mul(temp, temp, &C_pows[1], prec);
  acb_mul_si(temp, temp, 2, prec);
  acb_div_si(temp, temp, 3, prec);
  acb_add(&A_coefs[5], &A_coefs[5], temp, prec);

  /* Clear all power lists */
  _acb_vec_clear(A_pows, 2);
  _acb_vec_clear(B_pows, 3);
  _acb_vec_clear(C_pows, 2);

  /* A_coefs is now computed; compute C_coefs by hand */

  _acb_vec_set(C_coefs, A_coefs, 6);
  
  acb_pow_ui(temp, U, 2, prec);
  acb_mul(&C_coefs[0], &C_coefs[0], temp, prec);
  acb_pow_ui(temp, I10, 8, prec);
  acb_mul(&C_coefs[0], &C_coefs[0], temp, prec);
  
  acb_pow_ui(temp, I10, 10, prec);
  acb_mul(&C_coefs[1], &C_coefs[1], temp, prec);
  
  acb_pow_ui(temp, U, 8, prec);
  acb_mul(&C_coefs[2], &C_coefs[2], temp, prec);
  
  acb_pow_ui(temp, U, 4, prec);
  acb_mul(&C_coefs[3], &C_coefs[3], temp, prec);
  acb_pow_ui(temp, I10, 5, prec);
  acb_mul(&C_coefs[3], &C_coefs[3], temp, prec);
  
  acb_pow_ui(temp, U, 5, prec);
  acb_mul(&C_coefs[4], &C_coefs[4], temp, prec);
  acb_pow_ui(temp, I10, 4, prec);
  acb_mul(&C_coefs[4], &C_coefs[4], temp, prec);
  
  acb_pow_ui(temp, U, 1, prec);
  acb_mul(&C_coefs[5], &C_coefs[5], temp, prec);
  acb_pow_ui(temp, I10, 9, prec);
  acb_mul(&C_coefs[5], &C_coefs[5], temp, prec);
  
  _acb_vec_set(conic, C_coefs, 6);
  _acb_vec_clear(A_coefs, 6);
  _acb_vec_clear(C_coefs, 6);
  acb_clear(temp);
}
