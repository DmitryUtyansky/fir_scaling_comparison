/* Q15 instantiation of the FIR kernels: 16-bit samples, 32-bit products,
 * 64-bit wide accumulator. */
#include <math.h>
#include <stdint.h>

#include "fir_common.h"

#define FIR_PREFIX fir_q15
#define FIR_NAME   "Q15"
#define FIR_T      int16_t
#define FIR_TMIN   INT16_MIN
#define FIR_TMAX   INT16_MAX
#define FIR_T2     int32_t
#define FIR_T2U    uint32_t
#define FIR_T2MIN  INT32_MIN
#define FIR_T2MAX  INT32_MAX
#define FIR_T4     int64_t
#define FIR_W      16
#define FIR_FRAC   15

#include "fir_kernels.inc"
