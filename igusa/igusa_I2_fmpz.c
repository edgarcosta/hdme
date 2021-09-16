
#include "igusa.h"

void
igusa_I2_fmpz(fmpz_t res, const fmpz_t a0, const fmpz_t a1, const fmpz_t a2,
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
    a0_pows = _fmpz_vec_init(2);
    a1_pows = _fmpz_vec_init(2);
    a2_pows = _fmpz_vec_init(2);
    a3_pows = _fmpz_vec_init(3);
    a4_pows = _fmpz_vec_init(2);
    a5_pows = _fmpz_vec_init(2);
    a6_pows = _fmpz_vec_init(2);

/* Precompute all power lists */
    fmpz_one(&a0_pows[0]);
    fmpz_mul(&a0_pows[1], &a0_pows[0], a0);

    fmpz_one(&a1_pows[0]);
    fmpz_mul(&a1_pows[1], &a1_pows[0], a1);

    fmpz_one(&a2_pows[0]);
    fmpz_mul(&a2_pows[1], &a2_pows[0], a2);

    fmpz_one(&a3_pows[0]);
    fmpz_mul(&a3_pows[1], &a3_pows[0], a3);
    fmpz_mul(&a3_pows[2], &a3_pows[1], a3);

    fmpz_one(&a4_pows[0]);
    fmpz_mul(&a4_pows[1], &a4_pows[0], a4);

    fmpz_one(&a5_pows[0]);
    fmpz_mul(&a5_pows[1], &a5_pows[0], a5);

    fmpz_one(&a6_pows[0]);
    fmpz_mul(&a6_pows[1], &a6_pows[0], a6);

/* Add successive terms to res */
    fmpz_zero(res);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a3_pows[2]);
    fmpz_mul_si(temp, temp, 6);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a2_pows[1]);
    fmpz_mul(temp, temp, &a4_pows[1]);
    fmpz_mul_si(temp, temp, -16);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a1_pows[1]);
    fmpz_mul(temp, temp, &a5_pows[1]);
    fmpz_mul_si(temp, temp, 40);
    fmpz_add(res, res, temp);

    fmpz_one(temp);
    fmpz_mul(temp, temp, &a0_pows[1]);
    fmpz_mul(temp, temp, &a6_pows[1]);
    fmpz_mul_si(temp, temp, -240);
    fmpz_add(res, res, temp);

/* Clear all power lists */
    _fmpz_vec_clear(a0_pows, 2);
    _fmpz_vec_clear(a1_pows, 2);
    _fmpz_vec_clear(a2_pows, 2);
    _fmpz_vec_clear(a3_pows, 3);
    _fmpz_vec_clear(a4_pows, 2);
    _fmpz_vec_clear(a5_pows, 2);
    _fmpz_vec_clear(a6_pows, 2);

    fmpz_clear(temp);
}
