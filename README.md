# Fixed-point FIR scaling comparison

A stand-alone C project that implements a straightforward real FIR filter in
several fixed-point scaling styles and measures the accuracy of each of them
against a double-precision reference, in RMS SNR terms.

Everything is plain C11 plus `__int128` (a GCC extension used for the Q31
wide accumulator). No external dependencies beyond `libm`.

```
make            # build build/fir_compare
make test       # consistency self test
make run        # default summary tables (Q15 and Q31)
make sweep      # summary tables plus a scaling-shift sweep
```

## What is computed

For each output sample

```
y[n] = sum(k = 0 .. ntaps-1) h[k] * x[n + ntaps - 1 - k]
```

with `ntaps` = 256 and 1024 output samples by default. The input buffer holds
`len + ntaps - 1` samples, so every output is fully overlapped and there is no
start-up transient in the measurement.

Samples and coefficients are Qm numbers (m = 15 in a 16-bit container, or
m = 31 in a 32-bit container). The output has the same width as the input:
the accumulator is always scaled back to a single-width Qm result with
rounding and saturation.

### Scaling model

The exact filter output can reach the L1 gain of the filter, so the
fixed-point output is defined as

```
y_fx[n] ~= y_exact[n] / 2^G          G = ceil(log2(L1)),  L1 = sum |h|
```

`G` is the number of *guard bits*; it is derived from the actual quantised
coefficients and guarantees that the Qm output can never overflow. The
accumulator carries the product scale Q2m, so it has to be shifted right by

```
total = m + G
```

before it is written back as a Qm sample. Every variant applies exactly this
total shift, only split up differently between a shift in front of the
accumulator (`pre` in the report) and one after it (`post`), which always add
up to `total`. All nine results therefore live on the same scale and are
directly comparable.

The default `s` is `G`, the smallest shift that guarantees a 2x-wide
accumulator cannot overflow for full-scale input ("matched downscaling
shifts"). Both `s` and `G` can be forced from the command line.

### The nine variants

The listings below use these types and primitives, all of which exist in
`include/fir_kernels.inc` under longer names:

```c
T   x[], h[], y[];   /* sample width:  int16_t (Q15) / int32_t (Q31) */
T   acc;             /* single-wide accumulator                     */
T2  acc2;            /* 2x width:      int32_t      / int64_t        */
T4  acc4;            /* 4x width:      int64_t      / __int128       */

shift_rnd(v, s)  /* rounding right shift, in the type of v:
                    s ? (v + half(s)) >> s : v                        */
half(s)          /* rounding constant of a shift by s: s ? 1 << (s-1) : 0 */
sat(v)           /* clamp to the T range   (the output container)     */
sat2(v)          /* clamp to the T2 range  (the 2x accumulator)       */
wrap2(v)         /* two's complement wrap into T2, i.e. modulo 2^2W   */

xk               /* x[n + N-1 - k], the sample this tap feeds         */
r                /* the remaining shift, total - s                    */
```

**V1** — widening multiply into a 4x accumulator. Nothing is lost before the
final shift, so this is the accuracy upper bound.

```c
for (k = 0; k < N; k++)
    acc4 += (T2)xk * h[k];
y[n] = sat(shift_rnd(acc4, total));
```

**V2** — the input sample is pre-scaled before the multiply, otherwise V1.
The kernel hoists this shift into a single pass over the input buffer, which
is numerically identical and cheaper; V3 and V4 do the same.

```c
for (k = 0; k < N; k++)
    acc4 += (T2)shift_rnd(xk, s) * h[k];
y[n] = sat(shift_rnd(acc4, r));
```

**V3** — pre-scaled input, 2x accumulator, adds wrap around on overflow.

```c
for (k = 0; k < N; k++)
    acc2 = wrap2(acc2 + (T2)shift_rnd(xk, s) * h[k]);
y[n] = sat(wrap2(acc2 + half(r)) >> r);
```

**V4** — pre-scaled input, 2x accumulator, adds saturate on overflow.

```c
for (k = 0; k < N; k++)
    acc2 = sat2(acc2 + (T2)shift_rnd(xk, s) * h[k]);
y[n] = sat(sat2(acc2 + half(r)) >> r);
```

**V5** — full precision input, each product scaled down before a wrapping 2x
accumulator.

```c
for (k = 0; k < N; k++)
    acc2 = wrap2(acc2 + shift_rnd((T2)xk * h[k], s));
y[n] = sat(wrap2(acc2 + half(r)) >> r);
```

**V6** — as V5 with saturating adds.

```c
for (k = 0; k < N; k++)
    acc2 = sat2(acc2 + shift_rnd((T2)xk * h[k], s));
y[n] = sat(sat2(acc2 + half(r)) >> r);
```

**V7** — each product rounded straight onto the output LSB grid and summed at
the sample width. Nothing is left to shift, so `s` does not reach this one.

```c
for (k = 0; k < N; k++)
    acc = sat(acc + sat(shift_rnd((T2)xk * h[k], total)));
y[n] = acc;
```

**V8** — the same single-wide accumulator, but only `s` of the shift happens
in front of it.

```c
for (k = 0; k < N; k++)
    acc = sat(acc + sat(shift_rnd((T2)xk * h[k], s)));
y[n] = sat(shift_rnd(acc, r));
```

**V9** — the complementary split: the product is scaled by `r`, which at
matched shifts is exactly `m`, so it enters the accumulator as a plain Qm
sample.

```c
for (k = 0; k < N; k++)
    acc = sat(acc + sat(shift_rnd((T2)xk * h[k], r)));
y[n] = sat(shift_rnd(acc, s));
```

V7, V8 and V9 run the same kernel and differ only in where the total shift is
split around the accumulator. In V3 to V6 the rounding constant of the final
shift is written out because it has to travel through the same wrapping or
saturating adder as the rest of the sum; everywhere else `shift_rnd()` is
evaluated in a type that cannot overflow, so the plain form is exact.

### Reference and metrics

The reference is a double-precision FIR fed with the *quantised* input and
coefficients, converted to the output register domain by an exact power-of-two
division. The reported SNR therefore measures the accumulation and scaling
error of each variant, not the input quantisation error that all of them
share:

```
SNR = 10*log10( sum(ref^2) / sum((out - ref)^2) )
```

`rms err` and `max err` are in output LSBs. The harness also computes the
reference in `long double` and reports how far the double reference itself is
from exact ("double ref. floor"); for Q15 all products fit the 53-bit mantissa
exactly, so the reference is bit exact there.

Each variant is also compared against a closed-form error model (the `model`
column). Every variant ends with one output rounding, worth 1/12 LSB² of
noise variance; in addition

* pre-scaling the input by `s` bits injects 1/12 LSB² per input sample, which
  the filter amplifies by `sum(hq^2)` (its squared L2 gain in accumulator
  LSBs) before the final shift;
* scaling each product by `s` bits injects 1/12 LSB² per tap, i.e.
  `ntaps/12` at the accumulator, before the final shift. For V7 the products
  are rounded straight onto the output grid, so this term is `ntaps/12` at
  the output and there is no separate output rounding.

A variant fails if it is more than 3 dB below its model, which is what
`make test` and the exit status check. The status column reports
`OVERFLOW` / `SATURATE` / `CLIPPED` instead when the accumulator wrapped, the
accumulator clamped, or an output sample clamped: those are scaling failures
that the error model deliberately does not describe.

## Results for the default parameters

256 taps, 1024 samples, full-scale uniform white noise, random full-scale
coefficients (L1 = 128.2, so G = 8 guard bits), matched shifts `s = G = 8`:

### Q15 (16-bit samples, m = 15, total shift 23)

| # | variant | pre | post | SNR [dB] | model [dB] | rms err [LSB] | max err [LSB] |
|---|---|---|---|---|---|---|---|
| 1 | 4x acc, no pre-scale       | 0  | 23 | **67.4** | 67.5 | 0.290 | 0.50 |
| 2 | pre-scale + 4x acc         | 8  | 15 | 47.9 | 48.1 | 2.738 | 8.16 |
| 3 | pre-scale + 2x acc, wrap   | 8  | 15 | 47.9 | 48.1 | 2.738 | 8.16 |
| 4 | pre-scale + 2x acc, sat    | 8  | 15 | 47.9 | 48.1 | 2.738 | 8.16 |
| 5 | product >> s + 2x acc, wrap| 8  | 15 | **67.4** | 67.5 | 0.290 | 0.50 |
| 6 | product >> s + 2x acc, sat | 8  | 15 | **67.4** | 67.5 | 0.290 | 0.50 |
| 7 | product >> (s+r) + 1x acc  | 23 | 0  | 43.2 | 43.4 | 4.693 | 15.4 |
| 8 | product >> s + 1x acc      | 8  | 15 | 0.0 (`SATURATE`) | n/a | 680.7 | 2025 |
| 9 | product >> r + 1x acc      | 15 | 8  | 0.2 (`SATURATE`) | n/a | 663.9 | 2028 |

### Q31 (32-bit samples, m = 31, total shift 39)

| # | variant | pre | post | SNR [dB] | model [dB] | rms err [LSB] | max err [LSB] |
|---|---|---|---|---|---|---|---|
| 1 | 4x acc, no pre-scale       | 0  | 39 | **163.6** | 163.8 | 0.294 | 0.50 |
| 2 | pre-scale + 4x acc         | 8  | 31 | 144.4 | 144.4 | 2.700 | 8.61 |
| 3 | pre-scale + 2x acc, wrap   | 8  | 31 | 144.4 | 144.4 | 2.700 | 8.61 |
| 4 | pre-scale + 2x acc, sat    | 8  | 31 | 144.4 | 144.4 | 2.700 | 8.61 |
| 5 | product >> s + 2x acc, wrap| 8  | 31 | **163.6** | 163.8 | 0.294 | 0.50 |
| 6 | product >> s + 2x acc, sat | 8  | 31 | **163.6** | 163.8 | 0.294 | 0.50 |
| 7 | product >> (s+r) + 1x acc  | 39 | 0  | 139.7 | 139.7 | 4.604 | 15.9 |
| 8 | product >> s + 1x acc      | 8  | 31 | 0.0 (`SATURATE`) | n/a | 4.5e7 | 1.3e8 |
| 9 | product >> r + 1x acc      | 31 | 8  | 0.2 (`SATURATE`) | n/a | 4.4e7 | 1.3e8 |

Variants 1 to 7 see no accumulator overflow or saturation at `s = G` and no
output clipping; V2, V3 and V4 are bit exact with respect to each other, as
are V5 and V6. V8 and V9 clip, which is why no model value is shown for
them.

### Reading the table

* **The wide accumulator is worth ~19.5 dB over input pre-scaling.** V1 is
  limited only by the final output rounding (0.29 LSB rms, exactly
  `1/sqrt(12)`), and its 67.4 dB (Q15) is 33.7 dB below the 101.1 dB that the
  same rounding noise would give at full scale, simply because a random
  256-tap filter has an L1/L2 crest
  of 22.8 dB: the guard bits sized for the worst case leave the actual output
  at -33.7 dBFS. With the more realistic `--coefs lowpass` the same variant
  reaches 77.9 dB.
* **Scaling the products is far better than scaling the input, for the same
  overflow protection.** V5/V6 match the full-precision V1 exactly here, while
  V2/V3/V4 lose 19.5 dB. Pre-scaling discards `s` of the `m` significant bits
  of every input sample, and that error is then amplified by the L2 gain of
  the filter. Scaling the product by the same `s` bits removes `s` bits from a
  `2m`-bit product, which still leaves it well below the final output LSB
  (here the product-rounding noise lands 66 dB below the output rounding
  noise, so it is invisible).
* **With adequate guard bits, saturation costs nothing.** V4 == V3 and
  V6 == V5 bit exactly, because nothing ever overflows. The saturating logic
  only matters when the scaling is under-dimensioned.
* **Accumulating at single width costs 24 dB.** V7 rounds every product
  straight onto the output LSB grid, so instead of one output rounding there
  are `ntaps` of them, and the noise grows by `10*log10(ntaps) = 24.1 dB`
  (67.4 -> 43.2 dB, rms error 4.7 LSB). It never clips: the guard bits that
  bound the output also bound every partial sum, so a W-wide accumulator is
  sufficient here - it is accuracy, not range, that is lost.
* **A single-wide accumulator needs almost the whole downscale up front.**
  V8 keeps `total - s = 15` bits of the shift for after the loop, so it tries
  to store 22-bit products in a 16-bit accumulator: 381909 clamps and 0 dB
  SNR. V9 is the complementary split - it converts each product to a plain
  Q15 sample first, the textbook fractional multiply - and still fails, with
  50781 clamps, because the *sum* of 256 Q15 values needs the `G` guard bits
  that a 16-bit accumulator does not have. Only V7, which applies all 23 bits
  before the accumulator, survives at the default shift. Both degenerate to
  V7 when their split puts everything in front (V8 at `s = m + G`, V9 at
  `s = 0`), which the self test checks bit exactly.
* **Q31 behaves identically, 96 dB higher**, which is the expected 6 dB per
  extra bit.

### Under-dimensioned scaling: `--sweep`

`./build/fir_compare --type q15 --sweep` scans `s` and shows both failure
modes (`o` = the accumulator wrapped, `s` = it saturated):

```
    s  outsh      V2      V3      V4      V5      V6      V7      V8      V9
    0     23   67.4    -0.2o    0.8s   -0.2o    0.8s   43.2     0.0s   43.2
    1     22   67.4    -0.9o    3.3s   -0.9o    3.3s   43.2     0.0s   49.1
    2     21   67.3    -0.8o   10.9s   -0.8o   10.9s   43.2     0.0s   55.3
    3     20   67.1    67.1o   39.1s   67.4o   39.1s   43.2     0.0s   60.2
    4     19   66.2    66.2    66.2    67.4    67.4    43.2     0.0s   39.1s
    5     18   63.9    63.9    63.9    67.4    67.4    43.2     0.0s   10.9s
    6     17   59.7    59.7    59.7    67.4    67.4    43.2     0.0s    3.3s
    7     16   53.9    53.9    53.9    67.4    67.4    43.2     0.0s    0.8s
    8     15   47.9    47.9    47.9    67.4    67.4    43.2     0.0s    0.2s
    9     14   42.3    42.3    42.3    67.4    67.4    43.2     0.0s    0.1s
   ...
   16      7    0.0     0.0     0.0    67.4    67.4    43.2     0.8s    0.0s
   17      6    0.0     0.0     0.0    67.2    67.2    43.2     3.3s    0.0s
   18      5    0.0     0.0     0.0    66.4    66.4    43.2    10.9s    0.0s
   19      4    0.0     0.0     0.0    64.4    64.4    43.2    39.1s    0.0s
   20      3    0.0     0.0     0.0    60.2    60.2    43.2    60.2     0.0s
   21      2    0.0     0.0     0.0    55.3    55.3    43.2    55.3     0.0s
   22      1    0.0     0.0     0.0    49.1    49.1    43.2    49.1     0.0s
   23      0    0.0     0.0     0.0    43.2    43.2    43.2    43.2     0.0s
```

Several things are visible:

* Below `s = 4` the 2x accumulator can no longer hold the sum. The wrapping
  variants collapse to 0 dB (the error is as large as the signal), while the
  saturating ones degrade gradually.
* At `s = 3` the wrapping variants report wraps but still deliver full
  accuracy, whereas the saturating ones drop to 39 dB. Two's complement
  wrap-around is exact modulo 2^2W, so intermediate overflows cancel as long
  as the *final* sum fits; saturation destroys the value irreversibly. Once
  the final sum itself no longer fits (`s <= 2`), wrapping becomes the worse
  of the two.
* Above the minimum shift, V2/V3/V4 lose 6 dB per extra bit while V5/V6 stay
  flat until `s` grows so large that the product rounding reaches the output
  LSB (visible from `s = 17` on).
* V7 is a flat line at 43.2 dB: it always applies the full `m + G` shift to
  the products, so the sweep parameter does not reach it.
* V8 clips until `s` is within a few bits of `m + G`, then joins V5/V6 - the
  three are equal once the products fit the sample container, since a wider
  accumulator buys nothing when the values in it are already narrow. At the
  right edge (`s = m + G`, nothing left to shift) V5, V6, V7 and V8 all
  converge to the same 43.2 dB.
* V9 is the V8 column mirrored, since it shifts by `m + G - s` where V8 shifts
  by `s`. Its useful range is therefore at the *left* edge, and it shows what
  the single-wide accumulator can do at its best: leaving a few bits of the
  shift until after the loop lets the accumulator carry that many extra bits
  of the result, worth 6 dB each, until the running sum stops fitting.
  `s = 3` gives 60.2 dB, 17 dB better than V7 and only 7 dB below the wide
  accumulator, but it depends on the actual signal level rather than on a
  worst-case bound - one bit further (`s = 4`) and it clips to 39.1 dB. V7's
  split is the only one that is safe by construction.

## Sensitivity to the test data

Running the 3 signal kinds against the 3 coefficient sets (Q15, 256 taps,
default matched shifts) separates what depends on the data from what does
not. Measured SNR in dB:

| coefficients | L2 | signal | G | output | V1 | V2 | V5 | V7 | V9 | V1−V2 | V1−V7 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| random  | 9.28  | uniform | 8 | −33.7 | 67.4 | 47.9 | 67.4 | 43.2 | 0.2 `sat` | 19.5 | 24.2 |
| random  | 9.28  | gauss   | 8 | −40.7 | 60.4 | 40.8 | 60.4 | 36.3 | 1.1 `sat` | 19.5 | 24.1 |
| random  | 9.28  | chirp   | 8 | −31.5 | 69.5 | 50.1 | 69.5 | 45.2 | 0.2 `sat` | 19.4 | 24.3 |
| lowpass | 1.88  | uniform | 4 | −23.1 | 77.9 | 71.5 | 77.9 | 53.9 | 5.9 `sat` | 6.4 | 24.0 |
| lowpass | 1.88  | gauss   | 4 | −30.1 | 71.0 | 64.0 | 71.0 | 47.7 | 17.3 `sat` | 7.0 | 23.3 |
| lowpass | 1.88  | chirp   | 4 | −22.2 | 78.9 | 72.3 | 78.9 | 55.1 | 2.8 `sat` | 6.6 | 23.8 |
| l1norm  | 0.072 | uniform | 0 | −27.7 | 73.7 | 73.7 | 73.7 | 49.2 | 49.2 | 0.0 | 24.5 |
| l1norm  | 0.072 | gauss   | 0 | −34.7 | 66.0 | 66.0 | 66.0 | 42.6 | 42.6 | 0.0 | 23.5 |
| l1norm  | 0.072 | chirp   | 0 | −25.5 | 75.7 | 75.7 | 75.7 | 51.8 | 51.8 | 0.0 | 23.9 |

**The signal only moves the absolute level.** Within one coefficient set every
column shifts by the same amount, and that amount tracks the change in output
level to within a few tenths of a dB: gaussian noise at the default level sits
7.0 dB below uniform noise and costs 7.0 dB of SNR, the chirp sits 2.1 dB above
it and gains 2.1 dB. No ranking and no gap changes. This is expected, because
the noise of every variant is fixed in LSBs while only the signal power varies;
the residual tenths are the statistical spread of the realized rounding noise
over 1024 samples.

**The filter sets the pre-scaling penalty.** `V1 − V2` is `10*log10(1 + L2^2)`
at matched shifts, i.e. 19.40 / 6.58 / 0.02 dB for the three coefficient sets,
against 19.5 / 6.6 / 0.0 measured. The headline result that a wide accumulator
is worth ~19.5 dB is therefore a property of *random* coefficients: a real
low-pass filter only loses 6.6 dB, and a filter scaled to unit L1 gain loses
nothing at all.

**`V1 − V7` is 24 dB everywhere**, since `10*log10(1 + ntaps)` depends on
neither the signal nor the filter. The single-wide accumulator's penalty is
the most predictable number in the table.

**`l1norm` is a degenerate configuration.** With `L1 <= 1` there are no guard
bits, so the matched shift is `s = 0` and there is no scaling to compare:
variants 1 to 6 are bit identical, and V9 (whose `pre` becomes the full
`m + G`) is bit identical to V7. Only V8 still fails. Use it to check that
the harness reports "no scaling needed" correctly, not to compare variants.

**Overflow, unlike rounding, is strongly signal dependent.** The smallest
shift at which nothing wraps or clamps is well below the worst-case `G`, and
it moves with the signal:

| coefficients | G | uniform | gauss | chirp |
|---|---|---|---|---|
| random  | 8 | 4 | 2 | 4 |
| lowpass | 4 | 1 | 0 | 1 |
| l1norm  | 0 | 0 | 0 | 0 |

So does the answer to whether wrapping or saturation is better at the edge.
For random coefficients with uniform noise there is one shift (`s = 3`) where
the partial sums overflow but the totals still fit, and wrapping wins by 28 dB
(67.1 against 39.1). With the chirp, at the same shift, the coherent output is
large enough that the totals themselves overflow, and saturation wins instead
(21.1 against 4.3). With gaussian noise, whose rare peaks overflow before the
partial sums do, saturation wins at every shift where anything overflows at
all. The `wrap-is-harmless` result is a property of random-sign coefficients
driven by a non-peaky signal, not a general rule.

## Usage

```
usage: build/fir_compare [options]

  --taps N       number of FIR taps                     (default 256)
  --len N        number of output samples               (default 1024)
  --type T       q15 | q31 | both                       (default both)
  --shift S      scaling shift for variants 2..6        (default = guard bits)
  --guard G      output guard bits                      (default ceil(log2(L1)))
  --coefs C      random | lowpass | l1norm              (default random)
  --signal S     uniform | gauss | chirp                (default uniform)
  --level L      input peak level, linear 0..1          (default 1.0)
  --seed N       PRNG seed                              (default 1)
  --sweep        also sweep the scaling shift
  --selftest     run the consistency checks and exit
  --csv          emit the summary table as CSV
  --verbose      print extra counters
  --help         this text
```

The exit status is 0 when every variant meets its error model and 1 otherwise,
so both `make run` and `make test` are usable as regression checks.

## Self test

`make test` verifies, for both formats:

* V1 output is the correctly rounded exact result (peak error <= 0.5 LSB);
* V2 == V3 == V4 and V5 == V6 bit exactly when `s = G`;
* V8 == V7 bit exactly when `s = m + G`, and V9 == V7 when `s = 0`;
* every variant is within 3 dB of its error model (variants that clip or
  overflow are skipped, since the model does not describe them);
* with `s = 0` the 2x accumulator does report wrap-around, the saturating
  version does clamp, and clamping is the less damaging of the two;
* the impulse response reproduces the quantised coefficients.

## Layout

```
include/fir_common.h     types, variant list, public API
include/fir_kernels.inc  the kernels, instantiated once per data type
include/fir_ref.h        double / long double reference
include/testgen.h        reproducible signal and coefficient generation
include/metrics.h        SNR and error metrics
src/fir_q15.c            Q15 instantiation  (16/32/64-bit)
src/fir_q31.c            Q31 instantiation  (32/64/128-bit)
src/fir_common.c         variant names and descriptions
src/fir_ref.c            reference FIR
src/testgen.c            splitmix64 PRNG, signals, coefficient sets
src/metrics.c            SNR computation
src/main.c               harness, CLI, reporting, self test
```

The kernels are written once in `include/fir_kernels.inc` and instantiated
for each format by defining the sample type, the 2x and 4x wide types and the
number of fractional bits, so Q15 and Q31 run exactly the same code paths with
native types of the appropriate width.

Two implementation-defined behaviours are relied upon, both universal on GCC
and on any two's complement target: arithmetic right shift of negative
integers, and two's complement truncation when converting an out-of-range
unsigned value to a signed type. The latter is how the wrapping accumulator is
modelled without triggering signed-overflow undefined behaviour; overflow is
detected with the usual sign-bit test on the unsigned sum.
