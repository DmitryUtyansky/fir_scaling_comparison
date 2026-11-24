/*
 * fir_kernels.inc - fixed-point FIR kernels, instantiated once per data type.
 *
 * The including translation unit must define, before including this file:
 *
 *   FIR_PREFIX   symbol prefix, e.g. fir_q15
 *   FIR_NAME     format name string, e.g. "Q15"
 *   FIR_T        sample type                       (int16_t)
 *   FIR_TMIN     sample minimum                    (INT16_MIN)
 *   FIR_TMAX     sample maximum                    (INT16_MAX)
 *   FIR_T2       2x wide signed type               (int32_t)
 *   FIR_T2U      2x wide unsigned type             (uint32_t)
 *   FIR_T2MIN    2x minimum                        (INT32_MIN)
 *   FIR_T2MAX    2x maximum                        (INT32_MAX)
 *   FIR_T4       4x wide signed accumulator type   (int64_t)
 *   FIR_W        width of FIR_T in bits            (16)
 *   FIR_FRAC     fractional bits of FIR_T          (15)
 *
 * All names are #undef'd at the end so the file can be included repeatedly.
 *
 * Note: the code relies on two universally available implementation-defined
 * behaviours - arithmetic right shift of negative integers, and two's
 * complement truncation when converting an out-of-range unsigned value to a
 * signed type (used to model a wrapping accumulator without signed overflow).
 */

#define FIR_CAT2_(a, b) a##b
#define FIR_CAT_(a, b)  FIR_CAT2_(a, b)
#define FN(name)        FIR_CAT_(FIR_PREFIX, name)

#define FIR_SHIFT_SIGN (2 * FIR_W - 1)

/* ------------------------------------------------------------------ */
/* Primitives                                                          */
/* ------------------------------------------------------------------ */

/* Exact left shift, s >= 0, written as a multiply so that it stays defined
 * for negative v. */
static inline FIR_T4 FN(_shl4)(FIR_T4 v, int s)
{
    return v * ((FIR_T4)1 << s);
}

/*
 * Rounding constant added before a right shift of v by s > 0: half an output
 * LSB for round half up; for round half to even one less, plus the LSB of
 * the truncated result, which moves exact ties to the even neighbour and
 * leaves every other value alone; nothing for truncation (floor).  Never
 * larger than half, so it cannot overflow where half does not.
 */
static inline FIR_T2 FN(_rc2)(FIR_T2 v, int s, fir_rounding_t rnd)
{
    FIR_T2 half = (FIR_T2)1 << (s - 1);

    switch (rnd) {
    case FIR_ROUND_HALF_EVEN: return half - 1 + ((v >> s) & 1);
    case FIR_ROUND_NONE:      return 0;
    default:                  return half;
    }
}

static inline FIR_T4 FN(_rc4)(FIR_T4 v, int s, fir_rounding_t rnd)
{
    FIR_T4 half = (FIR_T4)1 << (s - 1);

    switch (rnd) {
    case FIR_ROUND_HALF_EVEN: return half - 1 + ((v >> s) & 1);
    case FIR_ROUND_NONE:      return 0;
    default:                  return half;
    }
}

/* Rounding right shift; a negative s is an exact left shift. */
static inline FIR_T4 FN(_rnd_shr4)(FIR_T4 v, int s, fir_rounding_t rnd)
{
    return (s > 0) ? ((v + FN(_rc4)(v, s, rnd)) >> s) : FN(_shl4)(v, -s);
}

/* Rounding right shift; s >= 0.  v is at most a full product,
 * |v| <= 2^(2W-2), so the rounding add only fits 2W for s <= 2W - 2. */
static inline FIR_T2 FN(_rnd_shr2)(FIR_T2 v, int s, fir_rounding_t rnd)
{
    if (s > 2 * FIR_W - 2)
        return (FIR_T2)FN(_rnd_shr4)((FIR_T4)v, s, rnd);
    return (s > 0) ? ((v + FN(_rc2)(v, s, rnd)) >> s) : v;
}

/* Clamp to the W-wide output container. */
static inline FIR_T FN(_sat_out)(FIR_T4 v, fir_stats_t *st)
{
    if (v > (FIR_T4)FIR_TMAX) { st->out_sat++; return FIR_TMAX; }
    if (v < (FIR_T4)FIR_TMIN) { st->out_sat++; return FIR_TMIN; }
    return (FIR_T)v;
}

/* Clamp to the 2W-wide output container. */
static inline FIR_T2 FN(_sat_out2)(FIR_T4 v, fir_stats_t *st)
{
    if (v > (FIR_T4)FIR_T2MAX) { st->out_sat++; return FIR_T2MAX; }
    if (v < (FIR_T4)FIR_T2MIN) { st->out_sat++; return FIR_T2MIN; }
    return (FIR_T2)v;
}

/* Write output sample n of a W-wide (wide == 0) or 2W-wide output buffer. */
static inline void FN(_store)(void *y, int n, int wide, FIR_T4 v,
                              fir_stats_t *st)
{
    if (wide)
        ((FIR_T2 *)y)[n] = FN(_sat_out2)(v, st);
    else
        ((FIR_T *)y)[n] = FN(_sat_out)(v, st);
}

/* 2W add that wraps around on overflow, counting the wraps. */
static inline FIR_T2 FN(_add_wrap2)(FIR_T2 a, FIR_T2 b, uint64_t *ovf)
{
    FIR_T2U ua = (FIR_T2U)a, ub = (FIR_T2U)b;
    FIR_T2U ur = (FIR_T2U)(ua + ub);

    *ovf += (uint64_t)((~(ua ^ ub) & (ua ^ ur)) >> FIR_SHIFT_SIGN);
    return (FIR_T2)ur;
}

/* 2W add that clamps on overflow, counting the clamps. */
static inline FIR_T2 FN(_add_sat2)(FIR_T2 a, FIR_T2 b, uint64_t *sat)
{
    FIR_T2U ua = (FIR_T2U)a, ub = (FIR_T2U)b;
    FIR_T2U ur = (FIR_T2U)(ua + ub);

    if ((~(ua ^ ub) & (ua ^ ur)) >> FIR_SHIFT_SIGN) {
        (*sat)++;
        return (a < 0) ? FIR_T2MIN : FIR_T2MAX;
    }
    return (FIR_T2)ur;
}

/* Clamp a scaled product into the W-wide sample container. */
static inline FIR_T FN(_sat1)(FIR_T2 v, uint64_t *sat)
{
    if (v > (FIR_T2)FIR_TMAX) { (*sat)++; return FIR_TMAX; }
    if (v < (FIR_T2)FIR_TMIN) { (*sat)++; return FIR_TMIN; }
    return (FIR_T)v;
}

/* W-wide add that clamps on overflow. */
static inline FIR_T FN(_add_sat1)(FIR_T a, FIR_T b, uint64_t *sat)
{
    return FN(_sat1)((FIR_T2)a + (FIR_T2)b, sat);
}

/*
 * Final scaling of a W-wide accumulator.  The rounding constant is added in
 * 2W, which is what a saturating W-wide rounding add followed by the shift
 * would produce once the result is clamped by FN(_store).  A negative s is
 * an exact left shift into a wider output.
 */
static inline FIR_T4 FN(_scale_out1)(FIR_T acc, int s, fir_rounding_t rnd)
{
    FIR_T2 v;

    if (s < 0)
        return FN(_shl4)((FIR_T4)acc, -s);
    v = (FIR_T2)acc + ((s > 0) ? FN(_rc2)((FIR_T2)acc, s, rnd) : 0);
    return (FIR_T4)(v >> s);
}

/* Final scaling of a 2W accumulator: rounding right shift (exact left shift
 * for a negative s); the caller saturates to the output. */
static inline FIR_T4 FN(_scale_out2)(FIR_T2 acc, int s, int saturating,
                                     fir_rounding_t rnd, fir_stats_t *st)
{
    FIR_T2 half;

    if (s < 0)
        return FN(_shl4)((FIR_T4)acc, -s);
    half = (s > 0) ? FN(_rc2)(acc, s, rnd) : 0;
    acc = saturating ? FN(_add_sat2)(acc, half, &st->acc_sat)
                     : FN(_add_wrap2)(acc, half, &st->acc_ovf);
    return (FIR_T4)(acc >> s);
}

/* ------------------------------------------------------------------ */
/* Kernels                                                             */
/* ------------------------------------------------------------------ */

/* Rounding right shift of a whole buffer (input pre-scaling). */
void FN(_prescale)(const FIR_T *src, FIR_T *dst, int n, int shift, fir_rounding_t rnd)
{
    int i;

    for (i = 0; i < n; i++)
        dst[i] = (FIR_T)FN(_rnd_shr2)((FIR_T2)src[i], shift, rnd);
}

/*
 * Coefficients of output n.  With h == NULL (squaring) every tap multiplies
 * its sample by the unscaled sample hx at the same position, so the
 * coefficients run backwards through hx exactly like the samples through x;
 * the coefficient of tap k is then hp[hs * k].
 */
static inline const FIR_T *FN(_coefs)(const FIR_T *h, const FIR_T *hx, int n,
                                      int ntaps, int *hs)
{
    *hs = h ? 1 : -1;
    return h ? h : hx + n + ntaps - 1;
}

/* V1 / V2: 2W widening products accumulated in a 4W accumulator. */
static void FN(_fir_acc4x)(const FIR_T *x, void *y, int wide, int nout,
                           const FIR_T *h, const FIR_T *hx, int ntaps,
                           int out_shift, fir_rounding_t rnd, fir_stats_t *st)
{
    int n, k, hs;

    for (n = 0; n < nout; n++) {
        const FIR_T *xp = x + n + ntaps - 1;
        const FIR_T *hp = FN(_coefs)(h, hx, n, ntaps, &hs);
        FIR_T4 acc = 0;

        for (k = 0; k < ntaps; k++)
            acc += (FIR_T4)((FIR_T2)xp[-k] * (FIR_T2)hp[hs * k]);

        FN(_store)(y, n, wide, FN(_rnd_shr4)(acc, out_shift, rnd), st);
    }
}

/* V3 / V4: 2W widening products accumulated in a 2W accumulator. */
static void FN(_fir_acc2x)(const FIR_T *x, void *y, int wide, int nout,
                           const FIR_T *h, const FIR_T *hx, int ntaps,
                           int out_shift, int saturating, fir_rounding_t rnd,
                           fir_stats_t *st)
{
    int n, k, hs;

    for (n = 0; n < nout; n++) {
        const FIR_T *xp = x + n + ntaps - 1;
        const FIR_T *hp = FN(_coefs)(h, hx, n, ntaps, &hs);
        FIR_T2 acc = 0;

        if (saturating) {
            for (k = 0; k < ntaps; k++)
                acc = FN(_add_sat2)(acc, (FIR_T2)xp[-k] * (FIR_T2)hp[hs * k],
                                    &st->acc_sat);
        } else {
            for (k = 0; k < ntaps; k++)
                acc = FN(_add_wrap2)(acc, (FIR_T2)xp[-k] * (FIR_T2)hp[hs * k],
                                     &st->acc_ovf);
        }
        FN(_store)(y, n, wide,
                   FN(_scale_out2)(acc, out_shift, saturating, rnd, st), st);
    }
}

/* V5 / V6: products scaled down (with rounding) before a 2W accumulator. */
static void FN(_fir_acc2x_prodshift)(const FIR_T *x, void *y, int wide,
                                     int nout, const FIR_T *h,
                                     const FIR_T *hx, int ntaps, int shift,
                                     int out_shift, int saturating, fir_rounding_t rnd,
                                     fir_stats_t *st)
{
    int n, k, hs;

    for (n = 0; n < nout; n++) {
        const FIR_T *xp = x + n + ntaps - 1;
        const FIR_T *hp = FN(_coefs)(h, hx, n, ntaps, &hs);
        FIR_T2 acc = 0;

        if (saturating) {
            for (k = 0; k < ntaps; k++) {
                FIR_T2 p = FN(_rnd_shr2)((FIR_T2)xp[-k] * (FIR_T2)hp[hs * k], shift,
                                         rnd);
                acc = FN(_add_sat2)(acc, p, &st->acc_sat);
            }
        } else {
            for (k = 0; k < ntaps; k++) {
                FIR_T2 p = FN(_rnd_shr2)((FIR_T2)xp[-k] * (FIR_T2)hp[hs * k], shift,
                                         rnd);
                acc = FN(_add_wrap2)(acc, p, &st->acc_ovf);
            }
        }
        FN(_store)(y, n, wide,
                   FN(_scale_out2)(acc, out_shift, saturating, rnd, st), st);
    }
}

/*
 * V7 / V8 / V9: products scaled down into the sample container and
 * accumulated at the sample width.  The caller decides how much of the total
 * shift happens before the accumulator (prod_shift) and how much after it
 * (post_shift); V7 puts all of it in front.
 */
static void FN(_fir_acc1x)(const FIR_T *x, void *y, int wide, int nout,
                           const FIR_T *h, const FIR_T *hx, int ntaps,
                           int prod_shift, int post_shift, fir_rounding_t rnd,
                           fir_stats_t *st)
{
    int n, k, hs;

    for (n = 0; n < nout; n++) {
        const FIR_T *xp = x + n + ntaps - 1;
        const FIR_T *hp = FN(_coefs)(h, hx, n, ntaps, &hs);
        FIR_T acc = 0;

        for (k = 0; k < ntaps; k++) {
            FIR_T2 p = FN(_rnd_shr2)((FIR_T2)xp[-k] * (FIR_T2)hp[hs * k],
                                     prod_shift, rnd);

            acc = FN(_add_sat1)(acc, FN(_sat1)(p, &st->acc_sat), &st->acc_sat);
        }
        FN(_store)(y, n, wide, FN(_scale_out1)(acc, post_shift, rnd), st);
    }
}

static void FN(_run_any)(fir_variant_t v, const FIR_T *x, FIR_T *scratch,
                         void *y, int wide, int nout, const FIR_T *h,
                         int ntaps, int shift, int out_shift,
                         fir_rounding_t rnd, fir_stats_t *st)
{
    int nx = nout + ntaps - 1;

    st->acc_ovf = st->acc_sat = st->out_sat = 0;

    switch (v) {
    case FIR_V1_ACC4X:
        FN(_fir_acc4x)(x, y, wide, nout, h, x, ntaps, out_shift, rnd, st);
        break;
    case FIR_V2_PRE_ACC4X:
        FN(_prescale)(x, scratch, nx, shift, rnd);
        FN(_fir_acc4x)(scratch, y, wide, nout, h, x, ntaps, out_shift, rnd, st);
        break;
    case FIR_V3_PRE_ACC2X_WRAP:
        FN(_prescale)(x, scratch, nx, shift, rnd);
        FN(_fir_acc2x)(scratch, y, wide, nout, h, x, ntaps, out_shift, 0, rnd, st);
        break;
    case FIR_V4_PRE_ACC2X_SAT:
        FN(_prescale)(x, scratch, nx, shift, rnd);
        FN(_fir_acc2x)(scratch, y, wide, nout, h, x, ntaps, out_shift, 1, rnd, st);
        break;
    case FIR_V5_POST_ACC2X_WRAP:
    case FIR_V10_PRODM_ACC2X_WRAP:
        FN(_fir_acc2x_prodshift)(x, y, wide, nout, h, x, ntaps, shift, out_shift, 0, rnd, st);
        break;
    case FIR_V6_POST_ACC2X_SAT:
        FN(_fir_acc2x_prodshift)(x, y, wide, nout, h, x, ntaps, shift, out_shift, 1, rnd, st);
        break;
    case FIR_V7_FULL_ACC1X_SAT:
    case FIR_V8_POST_ACC1X_SAT:
    case FIR_V9_PRE_ACC1X_SAT:
        /* The three differ only in how the caller splits the total shift. */
        FN(_fir_acc1x)(x, y, wide, nout, h, x, ntaps, shift, out_shift, rnd, st);
        break;
    default:
        break;
    }
}

void FN(_run)(fir_variant_t v, const FIR_T *x, FIR_T *scratch, FIR_T *y,
              int nout, const FIR_T *h, int ntaps, int shift, int out_shift,
              fir_rounding_t rnd, fir_stats_t *st)
{
    FN(_run_any)(v, x, scratch, y, 0, nout, h, ntaps, shift, out_shift, rnd,
                 st);
}

void FN(_run_wide)(fir_variant_t v, const FIR_T *x, FIR_T *scratch, FIR_T2 *y,
                   int nout, const FIR_T *h, int ntaps, int shift,
                   int out_shift, fir_rounding_t rnd, fir_stats_t *st)
{
    FN(_run_any)(v, x, scratch, y, 1, nout, h, ntaps, shift, out_shift, rnd,
                 st);
}

/* ------------------------------------------------------------------ */
/* Type-erased wrappers used by the harness                            */
/* ------------------------------------------------------------------ */

static void FN(_quantize)(const double *src, void *vdst, int n)
{
    const double scale = (double)((FIR_T4)1 << FIR_FRAC);
    FIR_T *dst = (FIR_T *)vdst;
    int i;

    for (i = 0; i < n; i++) {
        double v = src[i] * scale;
        double r = (v >= 0.0) ? floor(v + 0.5) : ceil(v - 0.5);

        if (r > (double)FIR_TMAX) r = (double)FIR_TMAX;
        if (r < (double)FIR_TMIN) r = (double)FIR_TMIN;
        dst[i] = (FIR_T)r;
    }
}

static void FN(_to_double)(const void *vsrc, double *dst, int n)
{
    const FIR_T *src = (const FIR_T *)vsrc;
    int i;

    for (i = 0; i < n; i++)
        dst[i] = (double)src[i];
}

static void FN(_wide_to_i64)(const void *vsrc, int64_t *dst, int n)
{
    const FIR_T2 *src = (const FIR_T2 *)vsrc;
    int i;

    for (i = 0; i < n; i++)
        dst[i] = (int64_t)src[i];
}

static void FN(_run_void)(fir_variant_t v, const void *x, void *scratch,
                          void *y, int nout, const void *h, int ntaps,
                          int shift, int out_shift, int out_wide,
                          fir_rounding_t rnd, fir_stats_t *st)
{
    FN(_run_any)(v, (const FIR_T *)x, (FIR_T *)scratch, y, out_wide, nout,
                 (const FIR_T *)h, ntaps, shift, out_shift, rnd, st);
}

const fir_type_ops_t FN(_ops) = {
    FIR_NAME,
    FIR_FRAC,
    FIR_W,
    sizeof(FIR_T),
    sizeof(FIR_T2),
    FN(_quantize),
    FN(_to_double),
    FN(_wide_to_i64),
    FN(_run_void)
};

#undef FN
#undef FIR_CAT_
#undef FIR_CAT2_
#undef FIR_SHIFT_SIGN
#undef FIR_PREFIX
#undef FIR_NAME
#undef FIR_T
#undef FIR_TMIN
#undef FIR_TMAX
#undef FIR_T2
#undef FIR_T2U
#undef FIR_T2MIN
#undef FIR_T2MAX
#undef FIR_T4
#undef FIR_W
#undef FIR_FRAC
