
void
igusa_R2_autogen(acb_t res, const acb_t I2, const acb_t I4, const acb_t I6,
                 const acb_t I10, slong prec)
{
    acb_t temp;
    acb_ptr I2_pows;
    acb_ptr I4_pows;
    acb_ptr I6_pows;
    acb_ptr I10_pows;
    acb_init(temp);

/* Init all power lists */
    I2_pows = _acb_vec_init(8);
    I4_pows = _acb_vec_init(8);
    I6_pows = _acb_vec_init(6);
    I10_pows = _acb_vec_init(4);

/* Precompute all power lists */
    acb_one(&I2_pows[0]);
    acb_mul(&I2_pows[1], &I2_pows[0], I2, prec);
    acb_mul(&I2_pows[2], &I2_pows[1], I2, prec);
    acb_mul(&I2_pows[3], &I2_pows[2], I2, prec);
    acb_mul(&I2_pows[4], &I2_pows[3], I2, prec);
    acb_mul(&I2_pows[5], &I2_pows[4], I2, prec);
    acb_mul(&I2_pows[6], &I2_pows[5], I2, prec);
    acb_mul(&I2_pows[7], &I2_pows[6], I2, prec);

    acb_one(&I4_pows[0]);
    acb_mul(&I4_pows[1], &I4_pows[0], I4, prec);
    acb_mul(&I4_pows[2], &I4_pows[1], I4, prec);
    acb_mul(&I4_pows[3], &I4_pows[2], I4, prec);
    acb_mul(&I4_pows[4], &I4_pows[3], I4, prec);
    acb_mul(&I4_pows[5], &I4_pows[4], I4, prec);
    acb_mul(&I4_pows[6], &I4_pows[5], I4, prec);
    acb_mul(&I4_pows[7], &I4_pows[6], I4, prec);

    acb_one(&I6_pows[0]);
    acb_mul(&I6_pows[1], &I6_pows[0], I6, prec);
    acb_mul(&I6_pows[2], &I6_pows[1], I6, prec);
    acb_mul(&I6_pows[3], &I6_pows[2], I6, prec);
    acb_mul(&I6_pows[4], &I6_pows[3], I6, prec);
    acb_mul(&I6_pows[5], &I6_pows[4], I6, prec);

    acb_one(&I10_pows[0]);
    acb_mul(&I10_pows[1], &I10_pows[0], I10, prec);
    acb_mul(&I10_pows[2], &I10_pows[1], I10, prec);
    acb_mul(&I10_pows[3], &I10_pows[2], I10, prec);

/* Add successive terms to res */
    acb_zero(res);

    acb_one(temp);
    acb_mul(temp, temp, &I2_pows[7], prec);
    acb_mul(temp, temp, &I4_pows[4], prec);
    acb_mul_si(temp, temp, 1, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &I2_pows[5], prec);
    acb_mul(temp, temp, &I4_pows[5], prec);
    acb_mul_si(temp, temp, 78, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &I2_pows[6], prec);
    acb_mul(temp, temp, &I4_pows[3], prec);
    acb_mul(temp, temp, &I6_pows[1], prec);
    acb_mul_si(temp, temp, -12, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &I2_pows[3], prec);
    acb_mul(temp, temp, &I4_pows[6], prec);
    acb_mul_si(temp, temp, -159, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &I2_pows[4], prec);
    acb_mul(temp, temp, &I4_pows[4], prec);
    acb_mul(temp, temp, &I6_pows[1], prec);
    acb_mul_si(temp, temp, -1332, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &I2_pows[5], prec);
    acb_mul(temp, temp, &I4_pows[2], prec);
    acb_mul(temp, temp, &I6_pows[2], prec);
    acb_mul_si(temp, temp, 54, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &I2_pows[6], prec);
    acb_mul(temp, temp, &I4_pows[2], prec);
    acb_mul(temp, temp, &I10_pows[1], prec);
    acb_mul_si(temp, temp, -972, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &I2_pows[1], prec);
    acb_mul(temp, temp, &I4_pows[7], prec);
    acb_mul_si(temp, temp, 80, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &I2_pows[2], prec);
    acb_mul(temp, temp, &I4_pows[5], prec);
    acb_mul(temp, temp, &I6_pows[1], prec);
    acb_mul_si(temp, temp, 1728, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &I2_pows[3], prec);
    acb_mul(temp, temp, &I4_pows[3], prec);
    acb_mul(temp, temp, &I6_pows[2], prec);
    acb_mul_si(temp, temp, 8910, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &I2_pows[4], prec);
    acb_mul(temp, temp, &I4_pows[1], prec);
    acb_mul(temp, temp, &I6_pows[3], prec);
    acb_mul_si(temp, temp, -108, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &I2_pows[4], prec);
    acb_mul(temp, temp, &I4_pows[3], prec);
    acb_mul(temp, temp, &I10_pows[1], prec);
    acb_mul_si(temp, temp, -77436, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &I2_pows[5], prec);
    acb_mul(temp, temp, &I4_pows[1], prec);
    acb_mul(temp, temp, &I6_pows[1], prec);
    acb_mul(temp, temp, &I10_pows[1], prec);
    acb_mul_si(temp, temp, 5832, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &I4_pows[6], prec);
    acb_mul(temp, temp, &I6_pows[1], prec);
    acb_mul_si(temp, temp, -384, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &I2_pows[1], prec);
    acb_mul(temp, temp, &I4_pows[4], prec);
    acb_mul(temp, temp, &I6_pows[2], prec);
    acb_mul_si(temp, temp, -6048, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &I2_pows[2], prec);
    acb_mul(temp, temp, &I4_pows[2], prec);
    acb_mul(temp, temp, &I6_pows[3], prec);
    acb_mul_si(temp, temp, -29376, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &I2_pows[3], prec);
    acb_mul(temp, temp, &I6_pows[4], prec);
    acb_mul_si(temp, temp, 81, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &I2_pows[2], prec);
    acb_mul(temp, temp, &I4_pows[4], prec);
    acb_mul(temp, temp, &I10_pows[1], prec);
    acb_mul_si(temp, temp, 592272, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &I2_pows[3], prec);
    acb_mul(temp, temp, &I4_pows[2], prec);
    acb_mul(temp, temp, &I6_pows[1], prec);
    acb_mul(temp, temp, &I10_pows[1], prec);
    acb_mul_si(temp, temp, 870912, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &I2_pows[4], prec);
    acb_mul(temp, temp, &I6_pows[2], prec);
    acb_mul(temp, temp, &I10_pows[1], prec);
    acb_mul_si(temp, temp, -8748, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &I2_pows[5], prec);
    acb_mul(temp, temp, &I10_pows[2], prec);
    acb_mul_si(temp, temp, 236196, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &I4_pows[3], prec);
    acb_mul(temp, temp, &I6_pows[3], prec);
    acb_mul_si(temp, temp, 6912, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &I2_pows[1], prec);
    acb_mul(temp, temp, &I4_pows[1], prec);
    acb_mul(temp, temp, &I6_pows[4], prec);
    acb_mul_si(temp, temp, 47952, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &I4_pows[5], prec);
    acb_mul(temp, temp, &I10_pows[1], prec);
    acb_mul_si(temp, temp, -41472, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &I2_pows[1], prec);
    acb_mul(temp, temp, &I4_pows[3], prec);
    acb_mul(temp, temp, &I6_pows[1], prec);
    acb_mul(temp, temp, &I10_pows[1], prec);
    acb_mul_si(temp, temp, -4743360, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &I2_pows[2], prec);
    acb_mul(temp, temp, &I4_pows[1], prec);
    acb_mul(temp, temp, &I6_pows[2], prec);
    acb_mul(temp, temp, &I10_pows[1], prec);
    acb_mul_si(temp, temp, -3090960, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &I2_pows[3], prec);
    acb_mul(temp, temp, &I4_pows[1], prec);
    acb_mul(temp, temp, &I10_pows[2], prec);
    acb_mul_si(temp, temp, 19245600, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &I6_pows[5], prec);
    acb_mul_si(temp, temp, -31104, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &I4_pows[2], prec);
    acb_mul(temp, temp, &I6_pows[2], prec);
    acb_mul(temp, temp, &I10_pows[1], prec);
    acb_mul_si(temp, temp, 9331200, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &I2_pows[1], prec);
    acb_mul(temp, temp, &I6_pows[3], prec);
    acb_mul(temp, temp, &I10_pows[1], prec);
    acb_mul_si(temp, temp, 3499200, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &I2_pows[1], prec);
    acb_mul(temp, temp, &I4_pows[2], prec);
    acb_mul(temp, temp, &I10_pows[2], prec);
    acb_mul_si(temp, temp, -507384000, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &I2_pows[2], prec);
    acb_mul(temp, temp, &I6_pows[1], prec);
    acb_mul(temp, temp, &I10_pows[2], prec);
    acb_mul_si(temp, temp, -104976000, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &I4_pows[1], prec);
    acb_mul(temp, temp, &I6_pows[1], prec);
    acb_mul(temp, temp, &I10_pows[2], prec);
    acb_mul_si(temp, temp, 2099520000, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

    acb_one(temp);
    acb_mul(temp, temp, &I10_pows[3], prec);
    acb_mul_si(temp, temp, 125971200000, prec);
    acb_div_si(temp, temp, 1, prec);
    acb_add(res, res, temp, prec);

/* Clear all power lists */
    _acb_vec_clear(I2_pows, 8);
    _acb_vec_clear(I4_pows, 8);
    _acb_vec_clear(I6_pows, 6);
    _acb_vec_clear(I10_pows, 4);

    acb_clear(temp);
}
