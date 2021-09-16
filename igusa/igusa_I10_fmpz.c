
#include "igusa.h"

void
igusa_I10_fmpz(fmpz_t res, const fmpz_t a0, const fmpz_t a1, const fmpz_t a2,
	       const fmpz_t a3, const fmpz_t a4, const fmpz_t a5,
	       const fmpz_t a6)
{
    fmpz_t temp;
    fmpz* a0_pows;
    fmpz* a1_pows;
    fmpz* a2_pows;
    fmpz* a3_pows;
    fmpz* a4_pows;
    fmpz* a5_pows;
    fmpz* a6_pows;
    fmpz_init(temp);

/* Init all power lists */
    a0_pows = _fmpz_vec_init(6);
    a1_pows = _fmpz_vec_init(7);
    a2_pows = _fmpz_vec_init(7);
    a3_pows = _fmpz_vec_init(7);
    a4_pows = _fmpz_vec_init(7);
    a5_pows = _fmpz_vec_init(7);
    a6_pows = _fmpz_vec_init(6);

/* Precompute all power lists */
    fmpz_one(&a0_pows[0]);
    fmpz_mul(&a0_pows[1], &a0_pows[0], a0);
    fmpz_mul(&a0_pows[2], &a0_pows[1], a0);
    fmpz_mul(&a0_pows[3], &a0_pows[2], a0);
    fmpz_mul(&a0_pows[4], &a0_pows[3], a0);
    fmpz_mul(&a0_pows[5], &a0_pows[4], a0);

    fmpz_one(&a1_pows[0]);
    fmpz_mul(&a1_pows[1], &a1_pows[0], a1);
    fmpz_mul(&a1_pows[2], &a1_pows[1], a1);
    fmpz_mul(&a1_pows[3], &a1_pows[2], a1);
    fmpz_mul(&a1_pows[4], &a1_pows[3], a1);
    fmpz_mul(&a1_pows[5], &a1_pows[4], a1);
    fmpz_mul(&a1_pows[6], &a1_pows[5], a1);

    fmpz_one(&a2_pows[0]);
    fmpz_mul(&a2_pows[1], &a2_pows[0], a2);
    fmpz_mul(&a2_pows[2], &a2_pows[1], a2);
    fmpz_mul(&a2_pows[3], &a2_pows[2], a2);
    fmpz_mul(&a2_pows[4], &a2_pows[3], a2);
    fmpz_mul(&a2_pows[5], &a2_pows[4], a2);
    fmpz_mul(&a2_pows[6], &a2_pows[5], a2);

    fmpz_one(&a3_pows[0]);
    fmpz_mul(&a3_pows[1], &a3_pows[0], a3);
    fmpz_mul(&a3_pows[2], &a3_pows[1], a3);
    fmpz_mul(&a3_pows[3], &a3_pows[2], a3);
    fmpz_mul(&a3_pows[4], &a3_pows[3], a3);
    fmpz_mul(&a3_pows[5], &a3_pows[4], a3);
    fmpz_mul(&a3_pows[6], &a3_pows[5], a3);

    fmpz_one(&a4_pows[0]);
    fmpz_mul(&a4_pows[1], &a4_pows[0], a4);
    fmpz_mul(&a4_pows[2], &a4_pows[1], a4);
    fmpz_mul(&a4_pows[3], &a4_pows[2], a4);
    fmpz_mul(&a4_pows[4], &a4_pows[3], a4);
    fmpz_mul(&a4_pows[5], &a4_pows[4], a4);
    fmpz_mul(&a4_pows[6], &a4_pows[5], a4);

    fmpz_one(&a5_pows[0]);
    fmpz_mul(&a5_pows[1], &a5_pows[0], a5);
    fmpz_mul(&a5_pows[2], &a5_pows[1], a5);
    fmpz_mul(&a5_pows[3], &a5_pows[2], a5);
    fmpz_mul(&a5_pows[4], &a5_pows[3], a5);
    fmpz_mul(&a5_pows[5], &a5_pows[4], a5);
    fmpz_mul(&a5_pows[6], &a5_pows[5], a5);

    fmpz_one(&a6_pows[0]);
    fmpz_mul(&a6_pows[1], &a6_pows[0], a6);
    fmpz_mul(&a6_pows[2], &a6_pows[1], a6);
    fmpz_mul(&a6_pows[3], &a6_pows[2], a6);
    fmpz_mul(&a6_pows[4], &a6_pows[3], a6);
    fmpz_mul(&a6_pows[5], &a6_pows[4], a6);

/* Add successive terms to res */
    fmpz_zero(res);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul_si(temp, temp, 1);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul_si(temp, temp, -4);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[3]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul_si(temp, temp, -4);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[3]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul_si(temp, temp, 18);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[4]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul_si(temp, temp, -27);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul_si(temp, temp, -4);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[4]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul_si(temp, temp, 16);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul_si(temp, temp, 18);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul_si(temp, temp, -80);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul_si(temp, temp, -6);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul_si(temp, temp, 144);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[4]);
    fmpz_mul(temp, temp, &a4_pows[4]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul_si(temp, temp, -27);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[4]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul_si(temp, temp, 144);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[4]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul_si(temp, temp, -128);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[4]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul_si(temp, temp, -192);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a4_pows[5]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul_si(temp, temp, 256);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[3]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul_si(temp, temp, -4);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[3]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul_si(temp, temp, 16);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[4]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul_si(temp, temp, 16);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[4]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul_si(temp, temp, -72);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[5]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul_si(temp, temp, 108);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul_si(temp, temp, 18);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[4]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul_si(temp, temp, -72);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul_si(temp, temp, -80);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul_si(temp, temp, 356);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[3]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul_si(temp, temp, 24);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[3]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul_si(temp, temp, -630);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul_si(temp, temp, -6);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul_si(temp, temp, 24);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[4]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul_si(temp, temp, 144);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul_si(temp, temp, -746);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul_si(temp, temp, 560);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul_si(temp, temp, 1020);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul_si(temp, temp, -36);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul_si(temp, temp, 160);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul_si(temp, temp, -1600);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[4]);
    fmpz_mul(temp, temp, &a5_pows[4]);
    fmpz_mul_si(temp, temp, -27);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[5]);
    fmpz_mul(temp, temp, &a5_pows[4]);
    fmpz_mul_si(temp, temp, 108);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[4]);
    fmpz_mul_si(temp, temp, 144);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[4]);
    fmpz_mul_si(temp, temp, -630);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[4]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[4]);
    fmpz_mul_si(temp, temp, -128);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[4]);
    fmpz_mul_si(temp, temp, 560);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[4]);
    fmpz_mul_si(temp, temp, 825);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[3]);
    fmpz_mul(temp, temp, &a5_pows[4]);
    fmpz_mul_si(temp, temp, -900);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[4]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[4]);
    fmpz_mul_si(temp, temp, -192);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[4]);
    fmpz_mul_si(temp, temp, 1020);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[4]);
    fmpz_mul_si(temp, temp, -900);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[4]);
    fmpz_mul_si(temp, temp, 160);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[4]);
    fmpz_mul_si(temp, temp, -2050);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[4]);
    fmpz_mul_si(temp, temp, 2250);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[4]);
    fmpz_mul_si(temp, temp, -50);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[4]);
    fmpz_mul_si(temp, temp, 2000);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[5]);
    fmpz_mul(temp, temp, &a5_pows[5]);
    fmpz_mul_si(temp, temp, 256);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[5]);
    fmpz_mul_si(temp, temp, -1600);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[5]);
    fmpz_mul_si(temp, temp, 2250);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[5]);
    fmpz_mul_si(temp, temp, 2000);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[5]);
    fmpz_mul_si(temp, temp, -3750);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[5]);
    fmpz_mul_si(temp, temp, -2500);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[4]);
    fmpz_mul(temp, temp, &a5_pows[6]);
    fmpz_mul_si(temp, temp, 3125);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -4);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 16);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[3]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 16);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[3]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -72);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[4]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 108);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a4_pows[4]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 16);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[4]);
    fmpz_mul(temp, temp, &a4_pows[4]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -64);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[4]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -72);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[4]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 320);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[4]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 24);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[4]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -576);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[4]);
    fmpz_mul(temp, temp, &a4_pows[5]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 108);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[5]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -576);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[5]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 512);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[5]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 768);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a4_pows[6]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -1024);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[3]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 18);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[3]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -72);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[4]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -72);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[4]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 324);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[5]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -486);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -80);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[4]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 320);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 356);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -1584);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[3]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -108);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[3]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 2808);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 24);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -96);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[4]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -630);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 3272);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -2496);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -4464);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a4_pows[4]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 144);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[4]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -640);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[4]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 6912);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -6);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[4]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 24);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[3]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 24);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[3]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -108);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[4]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 162);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[4]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 144);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[5]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -576);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -746);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 3272);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[4]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 560);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -2412);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -4536);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[3]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 3942);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[4]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 1020);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -5428);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 4816);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -682);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 10152);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -9720);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 248);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -10560);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -36);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[4]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 144);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[4]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 160);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -682);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -120);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -208);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 1980);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[3]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -1350);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[5]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -1600);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 9768);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -13040);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -12330);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 19800);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 15600);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[4]);
    fmpz_mul(temp, temp, &a5_pows[4]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 320);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[4]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -1700);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[4]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 1500);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[4]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 2250);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[4]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[4]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -22500);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[4]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -27);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[4]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 108);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[5]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 108);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[5]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -486);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[6]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 729);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 144);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[4]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -576);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[3]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -630);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[3]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 2808);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[4]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 162);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[4]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -4860);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[4]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -128);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[5]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 512);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 560);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -2496);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[4]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 825);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -4536);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 8208);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[3]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 5832);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[4]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -900);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 4816);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -4352);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -120);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -5760);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -8640);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[4]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -192);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[4]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 9216);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[4]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -192);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[5]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 768);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 1020);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -4464);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[4]);
    fmpz_mul(temp, temp, &a3_pows[3]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -900);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[3]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 3942);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[3]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 5832);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[4]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -6318);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 160);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[4]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -640);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[4]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -2050);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 10152);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -5760);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 1980);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -22896);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[3]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 21384);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[5]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 2250);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -13040);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 15264);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 16632);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -3456);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -21888);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[4]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -50);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 248);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[4]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -192);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[5]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 2000);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -12330);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 16632);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 15417);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -27540);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[4]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -1700);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 8748);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -6480);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -31320);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[4]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 43200);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 410);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -1800);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[4]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 27000);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[5]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, 256);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[6]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, -1024);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, -1600);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[4]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, 6912);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[4]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, 2250);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, -9720);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, -8640);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, -1350);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, 21384);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[4]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, -8748);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[4]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, 2000);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, -10560);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[4]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, 9216);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[5]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, -3750);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, 19800);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, -3456);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, -27540);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, 3888);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[4]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, 1500);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, -6480);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, -17280);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, 46656);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[4]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, -13824);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[5]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, -2500);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, 15600);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, -21888);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[4]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, 2250);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, -31320);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, 46656);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, 15552);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, -1800);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, 31968);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[4]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, -77760);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, 540);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[4]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, -32400);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[6]);
    fmpz_mul(temp, temp, &a6_pows[4]);
    fmpz_mul_si(temp, temp, 3125);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[4]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[4]);
    fmpz_mul_si(temp, temp, -22500);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[4]);
    fmpz_mul_si(temp, temp, 43200);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[4]);
    fmpz_mul_si(temp, temp, -13824);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[4]);
    fmpz_mul_si(temp, temp, 27000);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[4]);
    fmpz_mul_si(temp, temp, -77760);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[4]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[4]);
    fmpz_mul_si(temp, temp, 34992);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[4]);
    fmpz_mul_si(temp, temp, -32400);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[4]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[4]);
    fmpz_mul_si(temp, temp, 62208);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[4]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[4]);
    fmpz_mul_si(temp, temp, 38880);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[5]);
    fmpz_mul(temp, temp, &a6_pows[5]);
    fmpz_mul_si(temp, temp, -46656);
    fmpz_add(res, res, temp);

/* Clear all power lists */
    _fmpz_vec_clear(a0_pows, 6);
    _fmpz_vec_clear(a1_pows, 7);
    _fmpz_vec_clear(a2_pows, 7);
    _fmpz_vec_clear(a3_pows, 7);
    _fmpz_vec_clear(a4_pows, 7);
    _fmpz_vec_clear(a5_pows, 7);
    _fmpz_vec_clear(a6_pows, 6);

    fmpz_clear(temp);
}
