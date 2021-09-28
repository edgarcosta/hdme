
void
igusa_I4_autogen(acb_t res, const acb_t a0, const acb_t a1, const acb_t a2,
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
    a0_pows = _acb_vec_init(3);
    a1_pows = _acb_vec_init(3);
    a2_pows = _acb_vec_init(4);
    a3_pows = _acb_vec_init(3);
    a4_pows = _acb_vec_init(4);
    a5_pows = _acb_vec_init(3);
    a6_pows = _acb_vec_init(3);

/* Precompute all power lists */
    acb_one(&a0_pows[0]);
    acb_mul(&a0_pows[1], &a0_pows[0], a0, prec);
    acb_mul(&a0_pows[2], &a0_pows[1], a0, prec);

    acb_one(&a1_pows[0]);
    acb_mul(&a1_pows[1], &a1_pows[0], a1, prec);
    acb_mul(&a1_pows[2], &a1_pows[1], a1, prec);

    acb_one(&a2_pows[0]);
    acb_mul(&a2_pows[1], &a2_pows[0], a2, prec);
    acb_mul(&a2_pows[2], &a2_pows[1], a2, prec);
    acb_mul(&a2_pows[3], &a2_pows[2], a2, prec);

    acb_one(&a3_pows[0]);
    acb_mul(&a3_pows[1], &a3_pows[0], a3, prec);
    acb_mul(&a3_pows[2], &a3_pows[1], a3, prec);

    acb_one(&a4_pows[0]);
    acb_mul(&a4_pows[1], &a4_pows[0], a4, prec);
    acb_mul(&a4_pows[2], &a4_pows[1], a4, prec);
    acb_mul(&a4_pows[3], &a4_pows[2], a4, prec);

    acb_one(&a5_pows[0]);
    acb_mul(&a5_pows[1], &a5_pows[0], a5, prec);
    acb_mul(&a5_pows[2], &a5_pows[1], a5, prec);

    acb_one(&a6_pows[0]);
    acb_mul(&a6_pows[1], &a6_pows[0], a6, prec);
    acb_mul(&a6_pows[2], &a6_pows[1], a6, prec);

/* Add successive terms to res */
    acb_zero(res);

    acb_one(temp);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul_si(temp, temp, 4, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[2], prec);
    acb_mul_si(temp, temp, -12, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[3], prec);
    acb_mul_si(temp, temp, 48, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a2_pows[2], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul_si(temp, temp, -12, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul_si(temp, temp, 36, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul_si(temp, temp, 4, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul_si(temp, temp, -180, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul_si(temp, temp, -80, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[2], prec);
    acb_mul_si(temp, temp, 300, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a2_pows[3], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 48, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -180, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a3_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 324, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a1_pows[2], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, 300, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a2_pows[1], prec);
    acb_mul(temp, temp, &a4_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -504, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[1], prec);
    acb_mul(temp, temp, &a1_pows[1], prec);
    acb_mul(temp, temp, &a5_pows[1], prec);
    acb_mul(temp, temp, &a6_pows[1], prec);
    acb_mul_si(temp, temp, -540, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &a0_pows[2], prec);
    acb_mul(temp, temp, &a6_pows[2], prec);
    acb_mul_si(temp, temp, 1620, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

/* Clear all power lists */
    _acb_vec_clear(a0_pows, 3);
    _acb_vec_clear(a1_pows, 3);
    _acb_vec_clear(a2_pows, 4);
    _acb_vec_clear(a3_pows, 3);
    _acb_vec_clear(a4_pows, 4);
    _acb_vec_clear(a5_pows, 3);
    _acb_vec_clear(a6_pows, 3);

    acb_clear(temp);
}
