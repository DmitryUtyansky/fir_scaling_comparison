#include <math.h>

#include "metrics.h"

fir_metrics_t fir_metrics(const double *ref, const double *out, int n)
{
    fir_metrics_t m;
    double sig = 0.0, err = 0.0, peak = 0.0, maxe = 0.0;
    int i;

    for (i = 0; i < n; i++) {
        double e = out[i] - ref[i];
        double a = fabs(ref[i]);

        sig += ref[i] * ref[i];
        err += e * e;
        if (a > peak)
            peak = a;
        if (fabs(e) > maxe)
            maxe = fabs(e);
    }

    m.rms_sig  = (n > 0) ? sqrt(sig / (double)n) : 0.0;
    m.rms_err  = (n > 0) ? sqrt(err / (double)n) : 0.0;
    m.max_err  = maxe;
    m.peak_sig = peak;
    m.snr_db   = (err > 0.0) ? 10.0 * log10(sig / err) : INFINITY;
    return m;
}

double fir_ref_floor_db(const double *ref, const long double *ref_hi, int n)
{
    long double sig = 0.0L, err = 0.0L;
    int i;

    for (i = 0; i < n; i++) {
        long double e = (long double)ref[i] - ref_hi[i];

        sig += ref_hi[i] * ref_hi[i];
        err += e * e;
    }
    if (err <= 0.0L)
        return INFINITY;
    return 10.0 * log10((double)(sig / err));
}
