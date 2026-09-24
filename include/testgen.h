/*
 * testgen.h - reproducible test data generation.
 */
#ifndef TESTGEN_H
#define TESTGEN_H

#include <stdint.h>

typedef struct { uint64_t s; } rng_t;

void   rng_seed(rng_t *r, uint64_t seed);
double rng_uniform(rng_t *r);        /* uniform in [-1, 1)               */
double rng_gauss(rng_t *r);          /* zero mean, unit variance         */

typedef enum {
    SIG_UNIFORM = 0,   /* white uniform noise, full scale       */
    SIG_GAUSS,         /* white gaussian noise, -12 dBFS rms    */
    SIG_CHIRP          /* full scale linear sweep               */
} sig_kind_t;

typedef enum {
    COEF_RANDOM = 0,   /* white uniform coefficients, full scale */
    COEF_LOWPASS,      /* windowed-sinc low-pass, full scale     */
    COEF_L1NORM        /* random, normalised to unit L1 gain     */
} coef_kind_t;

int  sig_kind_parse(const char *s, sig_kind_t *out);
int  coef_kind_parse(const char *s, coef_kind_t *out);
const char *sig_kind_name(sig_kind_t k);
const char *coef_kind_name(coef_kind_t k);

/* Both generators produce values in (-1, 1). */
void gen_signal(rng_t *r, sig_kind_t kind, double level, double *x, int n);
void gen_coefs(rng_t *r, coef_kind_t kind, double *h, int ntaps);

#endif /* TESTGEN_H */
