
#include "igusa.h"

void
igusa_I6prime_fmpz(fmpz_t res, const fmpz_t a0, const fmpz_t a1,
		   const fmpz_t a2, const fmpz_t a3, const fmpz_t a4,
		   const fmpz_t a5, const fmpz_t a6)
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
    a0_pows = _fmpz_vec_init(4);
    a1_pows = _fmpz_vec_init(4);
    a2_pows = _fmpz_vec_init(5);
    a3_pows = _fmpz_vec_init(5);
    a4_pows = _fmpz_vec_init(5);
    a5_pows = _fmpz_vec_init(4);
    a6_pows = _fmpz_vec_init(4);

/* Precompute all power lists */
    fmpz_one(&a0_pows[0]);
    fmpz_mul(&a0_pows[1], &a0_pows[0], a0);
    fmpz_mul(&a0_pows[2], &a0_pows[1], a0);
    fmpz_mul(&a0_pows[3], &a0_pows[2], a0);

    fmpz_one(&a1_pows[0]);
    fmpz_mul(&a1_pows[1], &a1_pows[0], a1);
    fmpz_mul(&a1_pows[2], &a1_pows[1], a1);
    fmpz_mul(&a1_pows[3], &a1_pows[2], a1);

    fmpz_one(&a2_pows[0]);
    fmpz_mul(&a2_pows[1], &a2_pows[0], a2);
    fmpz_mul(&a2_pows[2], &a2_pows[1], a2);
    fmpz_mul(&a2_pows[3], &a2_pows[2], a2);
    fmpz_mul(&a2_pows[4], &a2_pows[3], a2);

    fmpz_one(&a3_pows[0]);
    fmpz_mul(&a3_pows[1], &a3_pows[0], a3);
    fmpz_mul(&a3_pows[2], &a3_pows[1], a3);
    fmpz_mul(&a3_pows[3], &a3_pows[2], a3);
    fmpz_mul(&a3_pows[4], &a3_pows[3], a3);

    fmpz_one(&a4_pows[0]);
    fmpz_mul(&a4_pows[1], &a4_pows[0], a4);
    fmpz_mul(&a4_pows[2], &a4_pows[1], a4);
    fmpz_mul(&a4_pows[3], &a4_pows[2], a4);
    fmpz_mul(&a4_pows[4], &a4_pows[3], a4);

    fmpz_one(&a5_pows[0]);
    fmpz_mul(&a5_pows[1], &a5_pows[0], a5);
    fmpz_mul(&a5_pows[2], &a5_pows[1], a5);
    fmpz_mul(&a5_pows[3], &a5_pows[2], a5);

    fmpz_one(&a6_pows[0]);
    fmpz_mul(&a6_pows[1], &a6_pows[0], a6);
    fmpz_mul(&a6_pows[2], &a6_pows[1], a6);
    fmpz_mul(&a6_pows[3], &a6_pows[2], a6);

/* Add successive terms to res */
    fmpz_zero(res);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul_si(temp, temp, 4);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul_si(temp, temp, -18);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul_si(temp, temp, 54);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[4]);
    fmpz_mul_si(temp, temp, 54);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[4]);
    fmpz_mul_si(temp, temp, -144);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul_si(temp, temp, -18);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul_si(temp, temp, 81);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[3]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul_si(temp, temp, -243);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul_si(temp, temp, 6);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul_si(temp, temp, -279);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul_si(temp, temp, 702);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul_si(temp, temp, 36);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a2_pows[4]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul_si(temp, temp, 54);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul_si(temp, temp, -279);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul_si(temp, temp, 216);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul_si(temp, temp, 405);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul_si(temp, temp, 624);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul_si(temp, temp, -1440);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul_si(temp, temp, -810);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul_si(temp, temp, 1350);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul_si(temp, temp, -1120);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul_si(temp, temp, 3600);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[3]);
    fmpz_mul_si(temp, temp, -3375);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 54);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -243);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[4]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 729);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a2_pows[4]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -144);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 702);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 405);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -3402);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -1440);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 2916);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 2754);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -5616);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 36);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -810);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 2754);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -2187);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 3600);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -11448);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 17010);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, 2160);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -8100);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 1350);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -5616);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[3]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -3375);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a3_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 17010);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -18954);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a1_pows[2]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, -8100);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 16524);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[2]);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[2]);
    fmpz_mul_si(temp, temp, 7290);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[3]);
    fmpz_mul(temp, temp, &a6_pows[3]);
    fmpz_mul_si(temp, temp, -14580);
    fmpz_add(res, res, temp);

/* Clear all power lists */
    _fmpz_vec_clear(a0_pows, 4);
    _fmpz_vec_clear(a1_pows, 4);
    _fmpz_vec_clear(a2_pows, 5);
    _fmpz_vec_clear(a3_pows, 5);
    _fmpz_vec_clear(a4_pows, 5);
    _fmpz_vec_clear(a5_pows, 4);
    _fmpz_vec_clear(a6_pows, 4);

    fmpz_clear(temp);
}
