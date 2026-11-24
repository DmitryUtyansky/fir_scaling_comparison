#include "fir_common.h"

static const struct {
    const char *name;
    const char *desc;
} variant_info[FIR_VARIANT_COUNT] = {
    { "acc4x",
      "2x widening mul, 4x accumulator, round+saturate to W" },
    { "pre + acc4x",
      "input >>s (rounded), 2x widening mul, 4x accumulator" },
    { "pre + acc2x wrap",
      "input >>s (rounded), 2x accumulator, wrapping adds" },
    { "pre + acc2x sat",
      "input >>s (rounded), 2x accumulator, saturating adds" },
    { "prod>>s + acc2x wrap",
      "2x mul, product >>s (rounded), 2x accumulator, wrapping adds" },
    { "prod>>s + acc2x sat",
      "2x mul, product >>s (rounded), 2x accumulator, saturating adds" },
    { "prod>>(s+r) + acc1x",
      "2x mul, product >>(s+r) to the output scale (rounded), 1x "
      "saturating accumulator, nothing left to shift" },
    { "prod>>s + acc1x",
      "2x mul, product >>s (rounded) clamped to W, 1x saturating "
      "accumulator, remaining shift >>r afterwards" },
    { "prod>>r + acc1x",
      "2x mul, product >>r (rounded) to a single-wide Qm value, 1x "
      "saturating accumulator, remaining shift >>s afterwards" },
    { "prod>>(inp_w-1)+acc2x wrap",
      "2x mul, product >>(W-1) (rounded) to a Qm value whatever the number "
      "of taps, 2x accumulator, wrapping adds, remaining shift afterwards" }
};

const char *fir_variant_name(fir_variant_t v)
{
    return (v < FIR_VARIANT_COUNT) ? variant_info[v].name : "?";
}

const char *fir_variant_desc(fir_variant_t v)
{
    return (v < FIR_VARIANT_COUNT) ? variant_info[v].desc : "?";
}

int fir_variant_uses_shift(fir_variant_t v)
{
    return v != FIR_V1_ACC4X;
}
