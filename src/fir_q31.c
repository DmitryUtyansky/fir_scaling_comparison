/* Q31 instantiation of the FIR kernels: 32-bit samples, 64-bit products,
 * 128-bit wide accumulator (GCC __int128). */
#include <math.h>
#include <stdint.h>

#include "fir_common.h"

#define FIR_PREFIX fir_q31
#define FIR_NAME   "Q31"
#define FIR_T      int32_t
#define FIR_TMIN   INT32_MIN
#define FIR_TMAX   INT32_MAX
#define FIR_T2     int64_t
#define FIR_T2U    uint64_t
#define FIR_T2MIN  INT64_MIN
#define FIR_T2MAX  INT64_MAX
#define FIR_T4     __int128
#define FIR_W      32
#define FIR_FRAC   31

#include "fir_kernels.inc"
