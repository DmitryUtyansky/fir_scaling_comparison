/*
 * fir_ref.h - floating-point reference FIR.
 */
#ifndef FIR_REF_H
#define FIR_REF_H

/*
 * y[n] = sum(k = 0 .. ntaps-1) h[k] * x[n + ntaps - 1 - k]
 *
 * x holds nout + ntaps - 1 samples, y holds nout samples.  Both x and h are
 * plain values; the harness feeds them as raw integer LSB counts so that the
 * result can be compared directly against the fixed-point output register.
 * h == NULL squares the input instead, as the kernels do.
 */
void fir_ref_f64(const double *x, const double *h, int ntaps, int nout,
                 double *y);

/* Same computation in long double, used to estimate the numerical noise
 * floor of the double-precision reference itself. */
void fir_ref_f80(const double *x, const double *h, int ntaps, int nout,
                 long double *y);

#endif /* FIR_REF_H */
