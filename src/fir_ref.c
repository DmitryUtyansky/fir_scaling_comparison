#include "fir_ref.h"

void fir_ref_f64(const double *x, const double *h, int ntaps, int nout,
                 double *y)
{
    int n, k;

    for (n = 0; n < nout; n++) {
        const double *xp = x + n + ntaps - 1;
        double acc = 0.0;

        for (k = 0; k < ntaps; k++)
            acc += h[k] * xp[-k];

        y[n] = acc;
    }
}

void fir_ref_f80(const double *x, const double *h, int ntaps, int nout,
                 long double *y)
{
    int n, k;

    for (n = 0; n < nout; n++) {
        const double *xp = x + n + ntaps - 1;
        long double acc = 0.0L;

        for (k = 0; k < ntaps; k++)
            acc += (long double)h[k] * (long double)xp[-k];

        y[n] = acc;
    }
}
