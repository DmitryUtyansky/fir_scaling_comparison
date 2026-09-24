#include <math.h>
#include <string.h>

#include "testgen.h"

/* splitmix64: small, fast, fully reproducible across platforms. */
static uint64_t rng_next(rng_t *r)
{
    uint64_t z = (r->s += 0x9E3779B97F4A7C15ULL);

    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

void rng_seed(rng_t *r, uint64_t seed)
{
    r->s = seed;
}

/* Uniform in [0, 1) with 53 significant bits. */
static double rng_unit(rng_t *r)
{
    return (double)(rng_next(r) >> 11) * (1.0 / 9007199254740992.0);
}

double rng_uniform(rng_t *r)
{
    return 2.0 * rng_unit(r) - 1.0;
}

double rng_gauss(rng_t *r)
{
    double u1 = rng_unit(r);
    double u2 = rng_unit(r);

    if (u1 < 1e-300)
        u1 = 1e-300;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

int sig_kind_parse(const char *s, sig_kind_t *out)
{
    if (!strcmp(s, "uniform")) { *out = SIG_UNIFORM; return 0; }
    if (!strcmp(s, "gauss"))   { *out = SIG_GAUSS;   return 0; }
    if (!strcmp(s, "chirp"))   { *out = SIG_CHIRP;   return 0; }
    return -1;
}

int coef_kind_parse(const char *s, coef_kind_t *out)
{
    if (!strcmp(s, "random"))  { *out = COEF_RANDOM;  return 0; }
    if (!strcmp(s, "lowpass")) { *out = COEF_LOWPASS; return 0; }
    if (!strcmp(s, "l1norm"))  { *out = COEF_L1NORM;  return 0; }
    return -1;
}

const char *sig_kind_name(sig_kind_t k)
{
    switch (k) {
    case SIG_UNIFORM: return "uniform white noise";
    case SIG_GAUSS:   return "gaussian white noise";
    case SIG_CHIRP:   return "linear chirp";
    }
    return "?";
}

const char *coef_kind_name(coef_kind_t k)
{
    switch (k) {
    case COEF_RANDOM:  return "random uniform, full scale";
    case COEF_LOWPASS: return "windowed-sinc low-pass, full scale";
    case COEF_L1NORM:  return "random uniform, unit L1 gain";
    }
    return "?";
}

void gen_signal(rng_t *r, sig_kind_t kind, double level, double *x, int n)
{
    int i;

    switch (kind) {
    case SIG_GAUSS:
        /* 4 sigma fits the full scale at the default level. */
        for (i = 0; i < n; i++) {
            double v = 0.25 * level * rng_gauss(r);

            if (v >= 1.0)  v = 0.999999;
            if (v <= -1.0) v = -0.999999;
            x[i] = v;
        }
        break;
    case SIG_CHIRP:
        for (i = 0; i < n; i++) {
            double t = (double)i / (double)(n > 1 ? n - 1 : 1);

            x[i] = level * 0.999 * sin(M_PI * (0.02 + 0.45 * t) * (double)i);
        }
        break;
    case SIG_UNIFORM:
    default:
        for (i = 0; i < n; i++)
            x[i] = level * rng_uniform(r);
        break;
    }
}

void gen_coefs(rng_t *r, coef_kind_t kind, double *h, int ntaps)
{
    int k;

    switch (kind) {
    case COEF_LOWPASS: {
        const double fc = 0.15;        /* cutoff, normalised to fs        */
        double peak = 0.0;

        for (k = 0; k < ntaps; k++) {
            double m = (double)k - 0.5 * (double)(ntaps - 1);
            double w = 0.54 - 0.46 * cos(2.0 * M_PI * (double)k /
                                         (double)(ntaps - 1));
            double s = (fabs(m) < 1e-12) ? 2.0 * fc
                                         : sin(2.0 * M_PI * fc * m) /
                                           (M_PI * m);
            h[k] = w * s;
            if (fabs(h[k]) > peak)
                peak = fabs(h[k]);
        }
        for (k = 0; k < ntaps; k++)      /* use the full coefficient range */
            h[k] *= 0.999 / peak;
        break;
    }
    case COEF_L1NORM: {
        double l1 = 0.0;

        for (k = 0; k < ntaps; k++) {
            h[k] = rng_uniform(r);
            l1 += fabs(h[k]);
        }
        if (l1 < 1e-30)
            l1 = 1.0;
        for (k = 0; k < ntaps; k++)
            h[k] *= 0.999 / l1;
        break;
    }
    case COEF_RANDOM:
    default:
        for (k = 0; k < ntaps; k++)
            h[k] = rng_uniform(r);
        break;
    }
}
