
#include "igusa.h"

void
igusa_I10_autogen(acb_t res, const acb_t a0, const acb_t a1, const acb_t a2,
                  const acb_t a3, const acb_t a4, const acb_t a5,
                  const acb_t a6, slong prec)
{
    acb_t temp;
    acb_ptr a0_pows;
    acb_ptr a1_pows;
    acb_ptr a2_pows;
    acb_ptr a3_pows;
    acb_ptr a4_pows;
    acb_ptr a5_pows;
    acb_ptr a6_pows;
    acb_init(temp);

/* Init all power lists */
    a0_pows = _acb_vec_init(6);
    a1_pows = _acb_vec_init(7);
    a2_pows = _acb_vec_init(7);
    a3_pows = _acb_vec_init(7);
    a4_pows = _acb_vec_init(7);
    a5_pows = _acb_vec_init(7);
    a6_pows = _acb_vec_init(6);

/* Precompute all power lists */
    acb_one(&a0_pows[0]);
    acb_mul(&a0_pows[1], &a0_pows[0], a0, prec);
    acb_mul(&a0_pows[2], &a0_pows[1], a0, prec);
    acb_mul(&a0_pows[3], &a0_pows[2], a0, prec);
    acb_mul(&a0_pows[4], &a0_pows[3], a0, prec);
    acb_mul(&a0_pows[5], &a0_pows[4], a0, prec);

    acb_one(&a1_pows[0]);
    acb_mul(&a1_pows[1], &a1_pows[0], a1, prec);
    acb_mul(&a1_pows[2], &a1_pows[1], a1, prec);
    acb_mul(&a1_pows[3], &a1_pows[2], a1, prec);
    acb_mul(&a1_pows[4], &a1_pows[3], a1, prec);
    acb_mul(&a1_pows[5], &a1_pows[4], a1, prec);
    acb_mul(&a1_pows[6], &a1_pows[5], a1, prec);

    acb_one(&a2_pows[0]);
    acb_mul(&a2_pows[1], &a2_pows[0], a2, prec);
    acb_mul(&a2_pows[2], &a2_pows[1], a2, prec);
    acb_mul(&a2_pows[3], &a2_pows[2], a2, prec);
    acb_mul(&a2_pows[4], &a2_pows[3], a2, prec);
    acb_mul(&a2_pows[5], &a2_pows[4], a2, prec);
    acb_mul(&a2_pows[6], &a2_pows[5], a2, prec);

    acb_one(&a3_pows[0]);
    acb_mul(&a3_pows[1], &a3_pows[0], a3, prec);
    acb_mul(&a3_pows[2], &a3_pows[1], a3, prec);
    acb_mul(&a3_pows[3], &a3_pows[2], a3, prec);
    acb_mul(&a3_pows[4], &a3_pows[3], a3, prec);
    acb_mul(&a3_pows[5], &a3_pows[4], a3, prec);
    acb_mul(&a3_pows[6], &a3_pows[5], a3, prec);

    acb_one(&a4_pows[0]);
    acb_mul(&a4_pows[1], &a4_pows[0], a4, prec);
    acb_mul(&a4_pows[2], &a4_pows[1], a4, prec);
    acb_mul(&a4_pows[3], &a4_pows[2], a4, prec);
    acb_mul(&a4_pows[4], &a4_pows[3], a4, prec);
    acb_mul(&a4_pows[5], &a4_pows[4], a4, prec);
    acb_mul(&a4_pows[6], &a4_pows[5], a4, prec);

    acb_one(&a5_pows[0]);
    acb_mul(&a5_pows[1], &a5_pows[0], a5, prec);
    acb_mul(&a5_pows[2], &a5_pows[1], a5, prec);
    acb_mul(&a5_pows[3], &a5_pows[2], a5, prec);
    acb_mul(&a5_pows[4], &a5_pows[3], a5, prec);
    acb_mul(&a5_pows[5], &a5_pows[4], a5, prec);
    acb_mul(&a5_pows[6], &a5_pows[5], a5, prec);

    acb_one(&a6_pows[0]);
    acb_mul(&a6_pows[1], &a6_pows[0], a6, prec);
    acb_mul(&a6_pows[2], &a6_pows[1], a6, prec);
    acb_mul(&a6_pows[3], &a6_pows[2], a6, prec);
    acb_mul(&a6_pows[4], &a6_pows[3], a6, prec);
    acb_mul(&a6_pows[5], &a6_pows[4], a6, prec);

/* Add successive terms to res */
    acb_zero(res);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul_si(temp, temp, -4, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[3], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul_si(temp, temp, -4, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[3], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul_si(temp, temp, 18, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[4], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul_si(temp, temp, -27, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[3], prec);
    acb_mul(temp, temp, &a4_pows[3], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul_si(temp, temp, -4, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[4], prec);
    acb_mul(temp, temp, &a4_pows[3], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul_si(temp, temp, 16, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[3], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul_si(temp, temp, 18, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[3], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul_si(temp, temp, -80, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[3], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul_si(temp, temp, -6, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[3], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul_si(temp, temp, 144, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[4], prec);
    acb_mul(temp, temp, &a4_pows[4], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul_si(temp, temp, -27, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[4], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul_si(temp, temp, 144, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[4], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul_si(temp, temp, -128, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[4], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul_si(temp, temp, -192, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a4_pows[5], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul_si(temp, temp, 256, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[3], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul_si(temp, temp, -4, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[3], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul_si(temp, temp, 16, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[4], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul_si(temp, temp, 16, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[4], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul_si(temp, temp, -72, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[5], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul_si(temp, temp, 108, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul_si(temp, temp, 18, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[4], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul_si(temp, temp, -72, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul_si(temp, temp, -80, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul_si(temp, temp, 356, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[3], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul_si(temp, temp, 24, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[3], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul_si(temp, temp, -630, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul_si(temp, temp, -6, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[3], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul_si(temp, temp, 24, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[4], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul_si(temp, temp, 144, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul_si(temp, temp, -746, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul_si(temp, temp, 560, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul_si(temp, temp, 1020, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a4_pows[3], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul_si(temp, temp, -36, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[3], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul_si(temp, temp, 160, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[3], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul_si(temp, temp, -1600, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[4], prec);
    acb_mul(temp, temp, &a5_pows[4], prec);
    acb_mul_si(temp, temp, -27, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[5], prec);
    acb_mul(temp, temp, &a5_pows[4], prec);
    acb_mul_si(temp, temp, 108, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[4], prec);
    acb_mul_si(temp, temp, 144, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[4], prec);
    acb_mul_si(temp, temp, -630, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[4], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[4], prec);
    acb_mul_si(temp, temp, -128, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[4], prec);
    acb_mul_si(temp, temp, 560, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[4], prec);
    acb_mul_si(temp, temp, 825, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[3], prec);
    acb_mul(temp, temp, &a5_pows[4], prec);
    acb_mul_si(temp, temp, -900, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[4], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[4], prec);
    acb_mul_si(temp, temp, -192, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[4], prec);
    acb_mul_si(temp, temp, 1020, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[3], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[4], prec);
    acb_mul_si(temp, temp, -900, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[4], prec);
    acb_mul_si(temp, temp, 160, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[4], prec);
    acb_mul_si(temp, temp, -2050, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[4], prec);
    acb_mul_si(temp, temp, 2250, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[4], prec);
    acb_mul_si(temp, temp, -50, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[4], prec);
    acb_mul_si(temp, temp, 2000, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[5], prec);
    acb_mul(temp, temp, &a5_pows[5], prec);
    acb_mul_si(temp, temp, 256, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[5], prec);
    acb_mul_si(temp, temp, -1600, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[5], prec);
    acb_mul_si(temp, temp, 2250, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[5], prec);
    acb_mul_si(temp, temp, 2000, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[5], prec);
    acb_mul_si(temp, temp, -3750, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[5], prec);
    acb_mul_si(temp, temp, -2500, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[4], prec);
    acb_mul(temp, temp, &a5_pows[6], prec);
    acb_mul_si(temp, temp, 3125, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[3], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -4, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[3], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 16, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[3], prec);
    acb_mul(temp, temp, &a4_pows[3], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 16, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[3], prec);
    acb_mul(temp, temp, &a4_pows[3], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -72, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[4], prec);
    acb_mul(temp, temp, &a4_pows[3], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 108, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[3], prec);
    acb_mul(temp, temp, &a4_pows[4], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 16, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[4], prec);
    acb_mul(temp, temp, &a4_pows[4], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -64, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[4], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -72, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[4], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 320, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[4], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 24, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[4], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -576, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[4], prec);
    acb_mul(temp, temp, &a4_pows[5], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 108, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[5], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -576, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[5], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 512, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[5], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 768, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a4_pows[6], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -1024, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[3], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 18, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[3], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -72, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[4], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -72, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[4], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 324, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[5], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -486, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -80, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[4], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 320, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 356, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -1584, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[3], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -108, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[3], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 2808, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[3], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 24, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[3], prec);
    acb_mul(temp, temp, &a4_pows[3], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -96, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[4], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[3], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -630, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[3], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 3272, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[3], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -2496, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[3], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -4464, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a4_pows[4], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 144, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[4], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -640, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[4], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 6912, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -6, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[4], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 24, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[3], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 24, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[3], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -108, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[4], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 162, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[4], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 144, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[5], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -576, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -746, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 3272, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[4], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 560, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -2412, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -4536, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[3], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 3942, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[4], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 1020, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -5428, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[3], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 4816, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -682, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 10152, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -9720, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[3], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 248, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[3], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -10560, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a2_pows[3], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -36, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[4], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 144, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[4], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 160, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -682, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -120, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -208, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 1980, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[3], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -1350, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[5], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -1600, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 9768, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -13040, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -12330, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 19800, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 15600, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[4], prec);
    acb_mul(temp, temp, &a5_pows[4], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 320, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[4], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -1700, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[4], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 1500, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[4], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 2250, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[4], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[4], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -22500, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[4], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, -27, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[4], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, 108, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[5], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, 108, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[5], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, -486, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[6], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, 729, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, 144, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[4], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, -576, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[3], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, -630, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[3], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, 2808, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[4], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, 162, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[4], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, -4860, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[4], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, -128, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[5], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, 512, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, 560, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, -2496, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[4], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, 825, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, -4536, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, 8208, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[3], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, 5832, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[4], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[3], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, -900, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[3], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, 4816, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[3], prec);
    acb_mul(temp, temp, &a4_pows[3], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, -4352, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[3], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, -120, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[3], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, -5760, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[3], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, -8640, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[4], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, -192, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[4], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, 9216, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[4], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, -192, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[5], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, 768, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, 1020, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, -4464, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[4], prec);
    acb_mul(temp, temp, &a3_pows[3], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, -900, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[3], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, 3942, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[3], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, 5832, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[4], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, -6318, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a2_pows[3], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, 160, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[4], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, -640, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[4], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, -2050, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, 10152, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, -5760, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, 1980, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, -22896, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[3], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, 21384, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[5], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, 2250, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, -13040, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, 15264, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, 16632, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, -3456, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[3], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, -21888, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[4], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, -50, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[3], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, 248, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[4], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, -192, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[5], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, 2000, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, -12330, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, 16632, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, 15417, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, -27540, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[4], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, -1700, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, 8748, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, -6480, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, -31320, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[4], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, 43200, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, 410, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, -1800, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[4], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[3], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, 27000, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[5], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, 256, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[6], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, -1024, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a2_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, -1600, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[4], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, 6912, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[4], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, 2250, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, -9720, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, -8640, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[3], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, -1350, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[3], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, 21384, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[4], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, -8748, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[4], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, 2000, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[3], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, -10560, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[4], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, 9216, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[5], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, -3750, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, 19800, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, -3456, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, -27540, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, 3888, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[4], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, 1500, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, -6480, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, -17280, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, 46656, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[4], prec);
    acb_mul(temp, temp, &a4_pows[3], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, -13824, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[5], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, -2500, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, 15600, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[3], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, -21888, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[4], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, 2250, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, -31320, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, 46656, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, 15552, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, -1800, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, 31968, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[4], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, -77760, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, 540, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[4], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[3], prec);
    acb_mul_si(temp, temp, -32400, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[6], prec);
    acb_mul(temp, temp, &a6_pows[4], prec);
    acb_mul_si(temp, temp, 3125, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[4], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[4], prec);
    acb_mul_si(temp, temp, -22500, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[4], prec);
    acb_mul_si(temp, temp, 43200, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a2_pows[3], prec);
    acb_mul(temp, temp, &a6_pows[4], prec);
    acb_mul_si(temp, temp, -13824, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a1_pows[3], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[4], prec);
    acb_mul_si(temp, temp, 27000, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[4], prec);
    acb_mul_si(temp, temp, -77760, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[4], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[4], prec);
    acb_mul_si(temp, temp, 34992, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[3], prec);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[4], prec);
    acb_mul_si(temp, temp, -32400, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[4], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[4], prec);
    acb_mul_si(temp, temp, 62208, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[4], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[4], prec);
    acb_mul_si(temp, temp, 38880, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[5], prec);
    acb_mul(temp, temp, &a6_pows[5], prec);
    acb_mul_si(temp, temp, -46656, prec);
    acb_add(res, res, temp, prec);

/* Clear all power lists */
    _acb_vec_clear(a0_pows, 6);
    _acb_vec_clear(a1_pows, 7);
    _acb_vec_clear(a2_pows, 7);
    _acb_vec_clear(a3_pows, 7);
    _acb_vec_clear(a4_pows, 7);
    _acb_vec_clear(a5_pows, 7);
    _acb_vec_clear(a6_pows, 6);

    acb_clear(temp);
}
