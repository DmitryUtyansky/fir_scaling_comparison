/*
 * fir_common.h - common definitions for the fixed-point FIR scaling comparison.
 */
#ifndef FIR_COMMON_H
#define FIR_COMMON_H

#include <stddef.h>
#include <stdint.h>

/*
 * Scaling / accumulation variants of the same real FIR.
 *
 * "W"  is the width of the input/output sample (16 for Q15, 32 for Q31),
 * "2W" is the width of the widening product, "4W" is twice that again.
 *
 *   V1  x*h -> 2W product, accumulated in a 4W accumulator, final rounding
 *       right shift + saturation back to W.  No information is ever lost
 *       before the final shift: this is the accuracy upper bound.
 *   V2  Input pre-scaled by a rounding right shift of `shift` bits, then
 *       identical to V1 (4W accumulator).  The final shift is reduced by
 *       `shift` so the output scale is unchanged.
 *   V3  Pre-scaled input, 2W accumulator, adds wrap around on overflow.
 *   V4  Pre-scaled input, 2W accumulator, adds saturate on overflow.
 *   V5  Full precision input, every 2W product is scaled down by a rounding
 *       right shift of `shift` bits before being added to a 2W accumulator,
 *       adds wrap around on overflow.
 *   V6  As V5, but the accumulator adds saturate.
 *   V7  Every 2W product is scaled down all the way to the output scale
 *       (`shift` plus the remaining shift, with rounding), clamped to W bits
 *       and accumulated in a W-wide saturating accumulator.  Nothing is left
 *       to shift after the loop.
 *   V8  As V7, but the products are only scaled down by `shift` before being
 *       clamped into the W-wide accumulator; the remaining shift is applied
 *       after the loop.
 *   V9  As V8 with the complementary split: the products are scaled down by
 *       the remaining shift `r` (which equals the fractional bit count at
 *       matched shifts, so the product becomes a plain Qm value) before the
 *       W-wide accumulator, and `shift` is applied after the loop.
 *   V10 As V5, but with the product shift fixed at W - 1 = m whatever the
 *       number of taps, so every product enters the accumulator as a Qm
 *       value; the remaining shift total - m is applied after the loop.
 *
 * V7, V8 and V9 are the same kernel; they only differ in how the total shift
 * is split around the accumulator, which is what the pre/post columns show.
 */
typedef enum {
    FIR_V1_ACC4X = 0,        /* 4x acc, no scaling before the accumulator   */
    FIR_V2_PRE_ACC4X,        /* pre-scaled input, 4x acc                    */
    FIR_V3_PRE_ACC2X_WRAP,   /* pre-scaled input, 2x acc, wrapping          */
    FIR_V4_PRE_ACC2X_SAT,    /* pre-scaled input, 2x acc, saturating        */
    FIR_V5_POST_ACC2X_WRAP,  /* scaled products,  2x acc, wrapping          */
    FIR_V6_POST_ACC2X_SAT,   /* scaled products,  2x acc, saturating        */
    FIR_V7_FULL_ACC1X_SAT,   /* fully scaled products, 1x acc, saturating   */
    FIR_V8_POST_ACC1X_SAT,   /* products >>s,  1x acc, saturating           */
    FIR_V9_PRE_ACC1X_SAT,    /* products >>r,  1x acc, saturating           */
    FIR_V10_PRODM_ACC2X_WRAP,/* products >>(W-1), 2x acc, wrapping          */
    FIR_VARIANT_COUNT
} fir_variant_t;

/* Event counters filled in by the kernels (diagnostics only). */
typedef struct {
    uint64_t acc_ovf;   /* wrap-around events in a 2W accumulator add */
    uint64_t acc_sat;   /* clamping events in a 2W accumulator add    */
    uint64_t out_sat;   /* clamping events when writing an output     */
} fir_stats_t;

/* Rounding of every right shift in the kernels: input pre-scaling, product
 * scaling and the final output shift. */
typedef enum {
    FIR_ROUND_HALF_UP = 0,   /* (v + half) >> s                         */
    FIR_ROUND_HALF_EVEN,     /* convergent: exact ties go to the even LSB */
    FIR_ROUND_NONE           /* truncation: v >> s, i.e. floor            */
} fir_rounding_t;

const char *fir_variant_name(fir_variant_t v);
const char *fir_variant_desc(fir_variant_t v);
int         fir_variant_uses_shift(fir_variant_t v);

/*
 * Typed kernel entry points.
 *
 * x        input, nout + ntaps - 1 samples (x[0] is the oldest sample)
 * scratch  buffer of nout + ntaps - 1 samples, used by the pre-scaling
 *          variants only, may be NULL otherwise
 * y        output, nout samples
 * h        coefficients, h[0] applies to the newest sample; NULL squares
 *          the input instead, i.e. every tap uses the (unscaled) sample it
 *          multiplies as its coefficient: y[n] = sum x[n + ntaps-1 - k]^2
 * shift    pre-scaling / product scaling shift (ignored by V1)
 * out_shift final rounding right shift applied to the accumulator; a
 *          negative value is an exact left shift (only useful with the
 *          2W-wide output of the _wide entry points)
 * rnd      rounding mode of all right shifts
 */
void fir_q15_run(fir_variant_t v, const int16_t *x, int16_t *scratch, int16_t *y,
                 int nout, const int16_t *h, int ntaps, int shift, int out_shift,
                 fir_rounding_t rnd, fir_stats_t *st);
void fir_q31_run(fir_variant_t v, const int32_t *x, int32_t *scratch, int32_t *y,
                 int nout, const int32_t *h, int ntaps, int shift, int out_shift,
                 fir_rounding_t rnd, fir_stats_t *st);

/* Same kernels with a 2W-wide output (Q31 for Q15 input, Q63 for Q31). */
void fir_q15_run_wide(fir_variant_t v, const int16_t *x, int16_t *scratch,
                      int32_t *y, int nout, const int16_t *h, int ntaps,
                      int shift, int out_shift, fir_rounding_t rnd,
                      fir_stats_t *st);
void fir_q31_run_wide(fir_variant_t v, const int32_t *x, int32_t *scratch,
                      int64_t *y, int nout, const int32_t *h, int ntaps,
                      int shift, int out_shift, fir_rounding_t rnd,
                      fir_stats_t *st);

/* Type-erased view of one fixed-point format, used by the test harness. */
typedef struct {
    const char *name;      /* "Q15" / "Q31"                                  */
    int         frac;      /* fractional bits (15 / 31)                      */
    int         width;     /* container width in bits (16 / 32)              */
    size_t      esize;     /* sizeof one sample                              */
    size_t      esize_wide;/* sizeof one 2W-wide output sample               */
    void      (*quantize)(const double *src, void *dst, int n);
    void      (*to_double)(const void *src, double *dst, int n); /* raw LSBs */
    void      (*wide_to_i64)(const void *src, int64_t *dst, int n); /* 2W out */
    void      (*run)(fir_variant_t v, const void *x, void *scratch, void *y,
                     int nout, const void *h, int ntaps, int shift,
                     int out_shift, int out_wide, fir_rounding_t rnd,
                     fir_stats_t *st);
} fir_type_ops_t;

extern const fir_type_ops_t fir_q15_ops;
extern const fir_type_ops_t fir_q31_ops;

#endif /* FIR_COMMON_H */
