/*
 * main.c - test harness for the fixed-point FIR scaling comparison.
 *
 * Scaling model
 * -------------
 * Input samples and coefficients are Qm numbers (m = 15 or 31) held in a
 * W-bit container (W = m + 1).  The exact filter output can be as large as
 * the L1 gain of the filter, so the fixed-point output is defined as
 *
 *     y_fx[n] ~= y_exact[n] / 2^G          G = ceil(log2(L1))   ("guard bits")
 *
 * which is guaranteed to stay inside the Qm range.  In LSB terms the
 * accumulator (which carries the product scale Q2m) therefore has to be
 * shifted right by
 *
 *     total = m + G
 *
 * before it is written back to a Qm output.  Every variant applies exactly
 * this total shift, split differently:
 *
 *     V1                shift = 0,  final shift = total
 *     V2, V3, V4        input >> s, final shift = total - s
 *     V5, V6            product >> s, final shift = total - s
 *
 * so all six results live on the same scale and can be compared directly.
 * The default s equals G, the smallest shift that guarantees a 2W
 * accumulator cannot overflow for full-scale input.
 *
 * The reference is a double-precision FIR fed with the *quantised* input and
 * coefficients, so the reported SNR measures the accumulation and scaling
 * error of each variant, not the input quantisation error that all of them
 * share.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fir_common.h"
#include "fir_ref.h"
#include "metrics.h"
#include "testgen.h"

#define PASS_MARGIN_DB 3.0   /* allowed shortfall against the model */

typedef struct {
    int taps;
    int len;
    int shift;
    int guard;
    int shift_set;
    int guard_set;
    uint64_t seed;
    double level;
    sig_kind_t sig;
    coef_kind_t coef;
    int do_q15;
    int do_q31;
    int sweep;
    int csv;
    int verbose;
} opts_t;

typedef struct {
    fir_metrics_t m;
    fir_stats_t   st;
    double        pred_db;
    int           shift;
    int           out_shift;
} result_t;

typedef struct {
    const fir_type_ops_t *ops;
    int taps, nout, nx;
    int guard, total;
    double l1, l2sq, sum_hq2;
    double ref_floor_db;
    double out_rms_dbfs, out_peak_dbfs;
    void   *xq, *hq, *scratch;
    void   *y[FIR_VARIANT_COUNT];
    double *ydbl;
    double *ref_reg;
} bench_t;

static void *xmalloc(size_t n)
{
    void *p = malloc(n ? n : 1);

    if (!p) {
        fprintf(stderr, "fir_compare: out of memory\n");
        exit(2);
    }
    return p;
}

static double db20(double x)
{
    return (x > 0.0) ? 20.0 * log10(x) : -INFINITY;
}

static int ceil_log2(double v)
{
    int g = 0;

    if (v <= 1.0)
        return 0;
    g = (int)ceil(log2(v) - 1e-12);
    return (g < 0) ? 0 : g;
}

/* ------------------------------------------------------------------ */
/* Benchmark setup                                                     */
/* ------------------------------------------------------------------ */

static void bench_init(bench_t *b, const fir_type_ops_t *ops, const opts_t *o)
{
    rng_t rng;
    double *xd, *hd, *hdq, *xdq, *ref_acc;
    long double *ref_hi;
    double scale, peak;
    int i, v;

    memset(b, 0, sizeof(*b));
    b->ops  = ops;
    b->taps = o->taps;
    b->nout = o->len;
    b->nx   = o->len + o->taps - 1;

    xd  = xmalloc((size_t)b->nx * sizeof(double));
    hd  = xmalloc((size_t)b->taps * sizeof(double));
    xdq = xmalloc((size_t)b->nx * sizeof(double));
    hdq = xmalloc((size_t)b->taps * sizeof(double));

    rng_seed(&rng, o->seed);
    gen_coefs(&rng, o->coef, hd, b->taps);
    gen_signal(&rng, o->sig, o->level, xd, b->nx);

    b->xq      = xmalloc((size_t)b->nx * ops->esize);
    b->hq      = xmalloc((size_t)b->taps * ops->esize);
    b->scratch = xmalloc((size_t)b->nx * ops->esize);
    for (v = 0; v < FIR_VARIANT_COUNT; v++)
        b->y[v] = xmalloc((size_t)b->nout * ops->esize);
    b->ydbl    = xmalloc((size_t)b->nout * sizeof(double));
    b->ref_reg = xmalloc((size_t)b->nout * sizeof(double));

    ops->quantize(xd, b->xq, b->nx);
    ops->quantize(hd, b->hq, b->taps);
    ops->to_double(b->xq, xdq, b->nx);
    ops->to_double(b->hq, hdq, b->taps);

    /* Gains of the quantised filter, in real (not LSB) units. */
    scale = ldexp(1.0, -ops->frac);
    for (i = 0; i < b->taps; i++) {
        double c = hdq[i] * scale;

        b->l1      += fabs(c);
        b->l2sq    += c * c;
        b->sum_hq2 += hdq[i] * hdq[i];
    }

    b->guard = o->guard_set ? o->guard : ceil_log2(b->l1);
    b->total = ops->frac + b->guard;

    /* Reference: exact sum of the integer products, scaled to the output
     * register domain.  The scaling is a power of two, hence exact. */
    ref_acc = xmalloc((size_t)b->nout * sizeof(double));
    ref_hi  = xmalloc((size_t)b->nout * sizeof(long double));
    fir_ref_f64(xdq, hdq, b->taps, b->nout, ref_acc);
    fir_ref_f80(xdq, hdq, b->taps, b->nout, ref_hi);

    peak = 0.0;
    for (i = 0; i < b->nout; i++) {
        b->ref_reg[i] = ldexp(ref_acc[i], -b->total);
        ref_hi[i]     = ldexpl(ref_hi[i], -b->total);
        if (fabs(b->ref_reg[i]) > peak)
            peak = fabs(b->ref_reg[i]);
    }
    b->ref_floor_db = fir_ref_floor_db(b->ref_reg, ref_hi, b->nout);

    {
        fir_metrics_t m = fir_metrics(b->ref_reg, b->ref_reg, b->nout);
        double full = ldexp(1.0, ops->frac);

        b->out_rms_dbfs  = db20(m.rms_sig / full);
        b->out_peak_dbfs = db20(peak / full);
    }

    free(xd); free(hd); free(xdq); free(hdq); free(ref_acc); free(ref_hi);
}

static void bench_free(bench_t *b)
{
    int v;

    free(b->xq); free(b->hq); free(b->scratch);
    for (v = 0; v < FIR_VARIANT_COUNT; v++)
        free(b->y[v]);
    free(b->ydbl); free(b->ref_reg);
}

/*
 * Predicted SNR.
 *
 * Every variant ends with one rounding of the output, contributing a
 * variance of 1/12 LSB^2.  In addition:
 *
 *  - pre-scaling the input by s bits injects an error of variance 1/12 LSB^2
 *    per input sample, which the filter amplifies by its squared L2 gain
 *    expressed in accumulator LSBs, sum(hq^2), before the final shift;
 *  - scaling every product by s bits injects 1/12 LSB^2 per tap, i.e.
 *    ntaps/12 at the accumulator, before the final shift.
 */
static double predict_snr_db(const bench_t *b, fir_variant_t v, int shift,
                             int out_shift)
{
    double sig = 0.0, noise, gain;
    int i;

    for (i = 0; i < b->nout; i++)
        sig += b->ref_reg[i] * b->ref_reg[i];
    sig /= (double)b->nout;

    /* The output rounding only costs something if there is a shift left. */
    noise = (out_shift > 0) ? 1.0 / 12.0 : 0.0;

    gain = ldexp(1.0, -2 * out_shift);
    if (shift > 0) {
        switch (v) {
        case FIR_V2_PRE_ACC4X:
        case FIR_V3_PRE_ACC2X_WRAP:
        case FIR_V4_PRE_ACC2X_SAT:
            noise += (1.0 / 12.0) * b->sum_hq2 * gain;
            break;
        case FIR_V5_POST_ACC2X_WRAP:
        case FIR_V6_POST_ACC2X_SAT:
        case FIR_V7_FULL_ACC1X_SAT:
        case FIR_V8_POST_ACC1X_SAT:
        case FIR_V9_PRE_ACC1X_SAT:
            noise += (1.0 / 12.0) * (double)b->taps * gain;
            break;
        default:
            break;
        }
    }
    return (noise > 0.0) ? 10.0 * log10(sig / noise) : INFINITY;
}

static void bench_run(bench_t *b, fir_variant_t v, int shift, result_t *r)
{
    int use_shift = fir_variant_uses_shift(v) ? shift : 0;

    /* V7 puts the whole shift in front of the accumulator, V9 uses the
     * complement of the swept shift; both keep pre + post = total. */
    if (v == FIR_V7_FULL_ACC1X_SAT)
        use_shift = b->total;
    else if (v == FIR_V9_PRE_ACC1X_SAT)
        use_shift = b->total - use_shift;

    r->shift     = use_shift;
    r->out_shift = b->total - use_shift;

    b->ops->run(v, b->xq, b->scratch, b->y[v], b->nout, b->hq, b->taps,
                use_shift, r->out_shift, &r->st);
    b->ops->to_double(b->y[v], b->ydbl, b->nout);
    r->m       = fir_metrics(b->ref_reg, b->ydbl, b->nout);
    r->pred_db = predict_snr_db(b, v, use_shift, r->out_shift);
}

/* ------------------------------------------------------------------ */
/* Reporting                                                           */
/* ------------------------------------------------------------------ */

static void fmt_db(char *buf, size_t n, double db)
{
    if (db == INFINITY)
        snprintf(buf, n, "%s", "exact");
    else
        snprintf(buf, n, "%.1f", db);
}

static const char *status_of(const result_t *r)
{
    if (r->st.acc_ovf)
        return "OVERFLOW";
    if (r->st.acc_sat)
        return "SATURATE";
    if (r->st.out_sat)
        return "CLIPPED";
    if (r->m.snr_db + PASS_MARGIN_DB < r->pred_db)
        return "FAIL";
    return "ok";
}

static int report(bench_t *b, const opts_t *o, result_t *res)
{
    int v, fails = 0;

    printf("\n");
    printf("================================================================"
           "================\n");
    printf("Fixed-point FIR scaling comparison - %s\n", b->ops->name);
    printf("================================================================"
           "================\n");
    printf("  format             : %s (%d-bit container, %d fractional bits)\n",
           b->ops->name, b->ops->width, b->ops->frac);
    printf("  output samples     : %d  (input buffer %d = %d + %d history)\n",
           b->nout, b->nx, b->nout, b->taps - 1);
    printf("  taps               : %d\n", b->taps);
    printf("  input signal       : %s, %.1f dBFS peak level\n",
           sig_kind_name(o->sig), db20(o->level));
    printf("  coefficients       : %s\n", coef_kind_name(o->coef));
    printf("  filter gain        : L1 = %.4f (%.2f dB), L2 = %.4f (%.2f dB)\n",
           b->l1, db20(b->l1), sqrt(b->l2sq), db20(sqrt(b->l2sq)));
    printf("  guard bits G       : %d%s  -> total right shift m+G = %d\n",
           b->guard, o->guard_set ? " (forced)" : " (= ceil(log2(L1)))",
           b->total);
    printf("  scaling shift s    : %d%s  (V7 always uses the full m+G, "
           "V9 the complement m+G-s)\n", res[FIR_V2_PRE_ACC4X].shift,
           o->shift_set ? " (forced)" : " (= G, matched downscaling)");
    printf("  reference output   : %.2f dBFS rms, %.2f dBFS peak"
           " (%.1f dB of the range is worst-case headroom)\n",
           b->out_rms_dbfs, b->out_peak_dbfs, -b->out_peak_dbfs);
    if (b->ref_floor_db == INFINITY)
        printf("  double ref. floor  : exact (all products fit the 53-bit "
               "mantissa)\n");
    else
        printf("  double ref. floor  : %.1f dB below signal\n",
               b->ref_floor_db);
    printf("----------------------------------------------------------------"
           "----------------\n");
    printf(" #  variant               pre   post     SNR    model   rms err"
           "  max err    ovf    sat  status\n");
    printf("                                        [dB]     [dB]     [LSB]"
           "    [LSB]\n");
    printf("----------------------------------------------------------------"
           "----------------\n");

    for (v = 0; v < FIR_VARIANT_COUNT; v++) {
        const result_t *r = &res[v];
        char snr[16], pred[16], sh[8];
        const char *st = status_of(r);

        fmt_db(snr, sizeof(snr), r->m.snr_db);
        if (strcmp(st, "ok") && strcmp(st, "FAIL"))
            snprintf(pred, sizeof(pred), "n/a");   /* model excludes clipping */
        else
            fmt_db(pred, sizeof(pred), r->pred_db);
        snprintf(sh, sizeof(sh), "%d", r->shift);

        printf(" %d  %-20s %3s %6d %7s %8s %9.4g %8.4g %6llu %6llu  %s\n",
               v + 1, fir_variant_name((fir_variant_t)v), sh, r->out_shift,
               snr, pred, r->m.rms_err, r->m.max_err,
               (unsigned long long)r->st.acc_ovf,
               (unsigned long long)r->st.acc_sat, st);

        if (!strcmp(st, "FAIL"))
            fails++;
    }
    printf("----------------------------------------------------------------"
           "----------------\n");
    for (v = 0; v < FIR_VARIANT_COUNT; v++)
        printf("  %d: %s\n", v + 1, fir_variant_desc((fir_variant_t)v));
    printf("  pre / post: right shift applied before / after the "
           "accumulator; pre + post = m + G = %d\n", b->total);
    printf("  status: ok = within %.0f dB of the error model; OVERFLOW / "
           "SATURATE / CLIPPED mean\n"
           "  the scaling is under-dimensioned, which the model does not "
           "describe.\n", PASS_MARGIN_DB);

    if (o->verbose) {
        printf("\n  output saturation events: ");
        for (v = 0; v < FIR_VARIANT_COUNT; v++)
            printf("V%d=%llu ", v + 1,
                   (unsigned long long)res[v].st.out_sat);
        printf("\n");
    }
    return fails;
}

static void report_csv(const bench_t *b, result_t *res, int header)
{
    int v;

    if (header)
        printf("type,taps,len,guard,variant,name,pre_shift,post_shift,snr_db,"
               "model_db,rms_err_lsb,max_err_lsb,acc_ovf,acc_sat,out_sat\n");
    for (v = 0; v < FIR_VARIANT_COUNT; v++) {
        const result_t *r = &res[v];

        printf("%s,%d,%d,%d,%d,%s,%d,%d,%.3f,%.3f,%.6f,%.3f,%llu,%llu,%llu\n",
               b->ops->name, b->taps, b->nout, b->guard, v + 1,
               fir_variant_name((fir_variant_t)v), r->shift, r->out_shift,
               r->m.snr_db, r->pred_db, r->m.rms_err, r->m.max_err,
               (unsigned long long)r->st.acc_ovf,
               (unsigned long long)r->st.acc_sat,
               (unsigned long long)r->st.out_sat);
    }
}

/* Sweep of the scaling shift, showing the overflow / precision trade-off. */
static void report_sweep(bench_t *b)
{
    int s, v;
    int smax = b->total < 24 ? b->total : 24;

    printf("\n  scaling shift sweep, SNR in dB "
           "('o' = accumulator wrapped, 's' = accumulator saturated)\n");
    printf("  ------------------------------------------------------------"
           "----\n");
    printf("    s  outsh");
    for (v = 1; v < FIR_VARIANT_COUNT; v++)
        printf("      V%d", v + 1);
    printf("\n");

    for (s = 0; s <= smax; s++) {
        printf("  %3d %6d", s, b->total - s);
        for (v = 1; v < FIR_VARIANT_COUNT; v++) {
            result_t r;
            char snr[16];
            char flag = ' ';

            bench_run(b, (fir_variant_t)v, s, &r);
            fmt_db(snr, sizeof(snr), r.m.snr_db);
            if (r.st.acc_ovf)
                flag = 'o';
            else if (r.st.acc_sat)
                flag = 's';
            printf(" %6s%c", snr, flag);
        }
        printf("\n");
    }
    printf("  ------------------------------------------------------------"
           "----\n");
    printf("  note: wrap-around is exact modulo 2^%d, so intermediate wraps\n"
           "        are harmless as long as the final sum fits - unlike\n"
           "        saturation, which destroys the result irreversibly.\n"
           "  note: V7 always shifts by m+G before the accumulator and V9 by\n"
           "        m+G-s, so the V9 column is the V8 column mirrored.\n",
           2 * b->ops->width);
}

/* ------------------------------------------------------------------ */
/* Self test                                                           */
/* ------------------------------------------------------------------ */

static int expect(int cond, const char *what)
{
    printf("  [%s] %s\n", cond ? "PASS" : "FAIL", what);
    return cond ? 0 : 1;
}

static int selftest(const fir_type_ops_t *ops, const opts_t *base)
{
    opts_t o = *base;
    bench_t b;
    result_t r[FIR_VARIANT_COUNT];
    int v, fails = 0;
    size_t ybytes;

    printf("\nself test - %s\n", ops->name);
    bench_init(&b, ops, &o);
    ybytes = (size_t)b.nout * ops->esize;

    for (v = 0; v < FIR_VARIANT_COUNT; v++)
        bench_run(&b, (fir_variant_t)v, o.shift_set ? o.shift : b.guard, &r[v]);

    fails += expect(r[FIR_V1_ACC4X].m.max_err <= 0.5 + 1e-9,
                    "V1 output is the correctly rounded exact result");
    fails += expect(r[FIR_V1_ACC4X].st.acc_ovf == 0 &&
                    r[FIR_V1_ACC4X].st.acc_sat == 0,
                    "V1 never overflows the 4x accumulator");
    fails += expect(memcmp(b.y[FIR_V2_PRE_ACC4X], b.y[FIR_V3_PRE_ACC2X_WRAP],
                           ybytes) == 0,
                    "V2 == V3 bit exact (2x accumulator does not overflow)");
    fails += expect(memcmp(b.y[FIR_V3_PRE_ACC2X_WRAP], b.y[FIR_V4_PRE_ACC2X_SAT],
                           ybytes) == 0,
                    "V3 == V4 bit exact (no saturation with s = G)");
    fails += expect(memcmp(b.y[FIR_V5_POST_ACC2X_WRAP], b.y[FIR_V6_POST_ACC2X_SAT],
                           ybytes) == 0,
                    "V5 == V6 bit exact (no saturation with s = G)");
    for (v = 0; v < FIR_VARIANT_COUNT; v++) {
        const char *st = status_of(&r[v]);
        char msg[128];

        if (strcmp(st, "ok") && strcmp(st, "FAIL")) {
            printf("  [skip] V%d SNR %.1f dB, %s: the model does not cover "
                   "under-dimensioned scaling\n", v + 1, r[v].m.snr_db, st);
            continue;
        }
        snprintf(msg, sizeof(msg), "V%d SNR %.1f dB within %.0f dB of the "
                 "model %.1f dB", v + 1, r[v].m.snr_db, PASS_MARGIN_DB,
                 r[v].pred_db);
        fails += expect(r[v].m.snr_db + PASS_MARGIN_DB >= r[v].pred_db, msg);
    }

    /* V8 and V9 both degenerate to V7 when their split puts the whole shift
     * in front of the accumulator. */
    {
        result_t rr;

        bench_run(&b, FIR_V7_FULL_ACC1X_SAT, b.total, &rr);
        bench_run(&b, FIR_V8_POST_ACC1X_SAT, b.total, &rr);
        fails += expect(memcmp(b.y[FIR_V7_FULL_ACC1X_SAT],
                               b.y[FIR_V8_POST_ACC1X_SAT], ybytes) == 0,
                        "V8 == V7 bit exact when s = m + G");
        bench_run(&b, FIR_V9_PRE_ACC1X_SAT, 0, &rr);
        fails += expect(memcmp(b.y[FIR_V7_FULL_ACC1X_SAT],
                               b.y[FIR_V9_PRE_ACC1X_SAT], ybytes) == 0,
                        "V9 == V7 bit exact when s = 0");
    }

    /* Deliberate under-scaling: the 2x accumulator must overflow, the
     * saturating variants must clamp instead of wrapping. */
    if (b.guard > 0) {
        result_t rw, rs;

        bench_run(&b, FIR_V3_PRE_ACC2X_WRAP, 0, &rw);
        bench_run(&b, FIR_V4_PRE_ACC2X_SAT, 0, &rs);
        fails += expect(rw.st.acc_ovf > 0,
                        "V3 detects wrap-around when s = 0 (no guard bits)");
        fails += expect(rs.st.acc_sat > 0,
                        "V4 saturates when s = 0 (no guard bits)");
        fails += expect(rs.m.snr_db > rw.m.snr_db,
                        "saturation degrades more gracefully than wrapping");
    }

    /* Impulse response: y must reproduce the coefficients. */
    {
        double *xd = xmalloc((size_t)b.nx * sizeof(double));
        double *yd = xmalloc((size_t)b.nout * sizeof(double));
        double *hdq = xmalloc((size_t)b.taps * sizeof(double));
        fir_stats_t st;
        double worst = 0.0;
        int i;

        for (i = 0; i < b.nx; i++)
            xd[i] = 0.0;
        xd[b.taps - 1] = ldexp(1.0, -1);              /* 0.5 full scale */
        ops->quantize(xd, b.xq, b.nx);
        ops->to_double(b.hq, hdq, b.taps);
        ops->run(FIR_V1_ACC4X, b.xq, b.scratch, b.y[0], b.nout, b.hq, b.taps,
                 0, b.total, &st);
        ops->to_double(b.y[0], yd, b.nout);
        for (i = 0; i < b.taps && i < b.nout; i++) {
            double want = ldexp(hdq[i] * ldexp(1.0, ops->frac - 1), -b.total);
            double e = fabs(yd[i] - want);

            if (e > worst)
                worst = e;
        }
        fails += expect(worst <= 0.5 + 1e-9,
                        "impulse response reproduces the coefficients");
        free(xd); free(yd); free(hdq);
    }

    bench_free(&b);
    return fails;
}

/* ------------------------------------------------------------------ */
/* CLI                                                                 */
/* ------------------------------------------------------------------ */

static void usage(const char *prog)
{
    printf(
"usage: %s [options]\n"
"\n"
"  --taps N       number of FIR taps                     (default 256)\n"
"  --len N        number of output samples               (default 1024)\n"
"  --type T       q15 | q31 | both                       (default both)\n"
"  --shift S      scaling shift for variants 2..6        (default = guard bits)\n"
"  --guard G      output guard bits                      (default ceil(log2(L1)))\n"
"  --coefs C      random | lowpass | l1norm              (default random)\n"
"  --signal S     uniform | gauss | chirp                (default uniform)\n"
"  --level L      input peak level, linear 0..1          (default 1.0)\n"
"  --seed N       PRNG seed                              (default 1)\n"
"  --sweep        also sweep the scaling shift\n"
"  --selftest     run the consistency checks and exit\n"
"  --csv          emit the summary table as CSV\n"
"  --verbose      print extra counters\n"
"  --help         this text\n", prog);
}

static int need_arg(int i, int argc, const char *opt)
{
    if (i >= argc) {
        fprintf(stderr, "fir_compare: %s requires an argument\n", opt);
        exit(2);
    }
    return 0;
}

int main(int argc, char **argv)
{
    opts_t o;
    int i, fails = 0, do_selftest = 0, csv_header_done = 0;

    memset(&o, 0, sizeof(o));
    o.taps   = 256;
    o.len    = 1024;
    o.seed   = 1;
    o.level  = 1.0;
    o.sig    = SIG_UNIFORM;
    o.coef   = COEF_RANDOM;
    o.do_q15 = 1;
    o.do_q31 = 1;

    for (i = 1; i < argc; i++) {
        const char *a = argv[i];

        if (!strcmp(a, "--taps"))       { need_arg(++i, argc, a); o.taps = atoi(argv[i]); }
        else if (!strcmp(a, "--len"))   { need_arg(++i, argc, a); o.len = atoi(argv[i]); }
        else if (!strcmp(a, "--shift")) { need_arg(++i, argc, a); o.shift = atoi(argv[i]); o.shift_set = 1; }
        else if (!strcmp(a, "--guard")) { need_arg(++i, argc, a); o.guard = atoi(argv[i]); o.guard_set = 1; }
        else if (!strcmp(a, "--seed"))  { need_arg(++i, argc, a); o.seed = strtoull(argv[i], NULL, 0); }
        else if (!strcmp(a, "--level")) { need_arg(++i, argc, a); o.level = atof(argv[i]); }
        else if (!strcmp(a, "--type")) {
            need_arg(++i, argc, a);
            o.do_q15 = !strcmp(argv[i], "q15") || !strcmp(argv[i], "both");
            o.do_q31 = !strcmp(argv[i], "q31") || !strcmp(argv[i], "both");
            if (!o.do_q15 && !o.do_q31) {
                fprintf(stderr, "fir_compare: unknown type '%s'\n", argv[i]);
                return 2;
            }
        } else if (!strcmp(a, "--coefs")) {
            need_arg(++i, argc, a);
            if (coef_kind_parse(argv[i], &o.coef)) {
                fprintf(stderr, "fir_compare: unknown coefficient set '%s'\n", argv[i]);
                return 2;
            }
        } else if (!strcmp(a, "--signal")) {
            need_arg(++i, argc, a);
            if (sig_kind_parse(argv[i], &o.sig)) {
                fprintf(stderr, "fir_compare: unknown signal '%s'\n", argv[i]);
                return 2;
            }
        }
        else if (!strcmp(a, "--sweep"))    o.sweep = 1;
        else if (!strcmp(a, "--csv"))      o.csv = 1;
        else if (!strcmp(a, "--verbose"))  o.verbose = 1;
        else if (!strcmp(a, "--selftest")) do_selftest = 1;
        else if (!strcmp(a, "--help") || !strcmp(a, "-h")) { usage(argv[0]); return 0; }
        else {
            fprintf(stderr, "fir_compare: unknown option '%s'\n", a);
            usage(argv[0]);
            return 2;
        }
    }

    if (o.taps < 1 || o.len < 1) {
        fprintf(stderr, "fir_compare: --taps and --len must be >= 1\n");
        return 2;
    }
    if (o.level <= 0.0 || o.level > 1.0) {
        fprintf(stderr, "fir_compare: --level must be in (0, 1]\n");
        return 2;
    }

    if (do_selftest) {
        if (o.do_q15)
            fails += selftest(&fir_q15_ops, &o);
        if (o.do_q31)
            fails += selftest(&fir_q31_ops, &o);
        printf("\nself test: %s (%d failure%s)\n", fails ? "FAILED" : "PASSED",
               fails, fails == 1 ? "" : "s");
        return fails ? 1 : 0;
    }

    for (i = 0; i < 2; i++) {
        const fir_type_ops_t *ops = i ? &fir_q31_ops : &fir_q15_ops;
        result_t res[FIR_VARIANT_COUNT];
        bench_t b;
        int shift, v;

        if (i == 0 && !o.do_q15) continue;
        if (i == 1 && !o.do_q31) continue;

        bench_init(&b, ops, &o);
        shift = o.shift_set ? o.shift : b.guard;
        if (shift < 0 || shift > b.total) {
            fprintf(stderr, "fir_compare: --shift must be in [0, %d] for %s\n",
                    b.total, ops->name);
            bench_free(&b);
            return 2;
        }
        for (v = 0; v < FIR_VARIANT_COUNT; v++)
            bench_run(&b, (fir_variant_t)v, shift, &res[v]);

        if (o.csv)
            report_csv(&b, res, !csv_header_done++);
        else
            fails += report(&b, &o, res);

        if (o.sweep && !o.csv)
            report_sweep(&b);

        bench_free(&b);
    }

    if (!o.csv)
        printf("\n%s\n", fails ? "RESULT: FAIL - a variant is worse than its "
                                "error model" :
                                "RESULT: PASS - every variant matches its "
                                "error model");
    return fails ? 1 : 0;
}
