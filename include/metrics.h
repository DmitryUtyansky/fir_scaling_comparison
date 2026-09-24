/*
 * metrics.h - accuracy metrics, all expressed in output LSBs.
 */
#ifndef METRICS_H
#define METRICS_H

typedef struct {
    double snr_db;       /* 10*log10(sum ref^2 / sum err^2)   */
    double rms_err;      /* rms error, output LSB             */
    double max_err;      /* peak absolute error, output LSB   */
    double rms_sig;      /* rms of the reference, output LSB  */
    double peak_sig;     /* peak of the reference, output LSB */
} fir_metrics_t;

/* ref and out are both in the output register domain (LSB counts). */
fir_metrics_t fir_metrics(const double *ref, const double *out, int n);

/* Numerical noise floor of the double reference, in dB below the signal. */
double fir_ref_floor_db(const double *ref, const long double *ref_hi, int n);

#endif /* METRICS_H */
