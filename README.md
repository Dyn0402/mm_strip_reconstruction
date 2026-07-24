# mm_strip_reconstruction
Reconstruction algorithms build specifically for Micromegas 2D strip detectors built at CEA Saclay.

## Build

Use the **Release** build for processing — it is ~12× faster than the debug
build (which was unknowingly used for all processing before 2026-07-24):

```
cmake -S . -B cmake-build-release -DCMAKE_BUILD_TYPE=Release
cmake --build cmake-build-release
```

`orchestrator/process_run.py` points at `cmake-build-release/`. Reference
timing: one 13k-event non-ZS FEU file (32 samples × 512 ch) analyzes in ~12 s
(was ~160 s debug); remaining time is ~⅓ ROOT decompression, ~⅓ CNS medians.

## Pipeline

`decode` (raw fdf → `nt` sample tree) → `analyze_waveforms` (→ per-FEU `hits`
tree) → `combine_feus_hits` (→ combined hits with `feu` branch) →
`clusterize1d`. `orchestrator/process_run.py` drives the chain for a run.

## Waveform analysis (`waveform_analysis/`)

```
analyze_waveforms <input.root> [output.root] [pedestal.root] [--tps <ns>] [--cns <0|1>]
                  [--thr <sigma>] [--mf <samples>]
```

Per event: pedestal subtraction → densify zero-suppressed samples → optional
common-noise subtraction (median per 64-channel block per sample; ON by
default) → pulse finding → one `hits` entry per pulse.

Pedestals: mean = raw baseline per channel; RMS = noise measured **after** CNS
(the data path is CNS-subtracted, so thresholds must be calibrated on post-CNS
noise). Without a pedestal file the data is assumed zero-suppressed
(baseline 256, RMS 1). CNS on zero-suppressed (sparse) data is signal-biased
(the block median *is* signal → ~8× hit loss), and ZS data is already
common-mode corrected on the FEU, so the analyzer **auto-detects ZS input**
(median distinct channels/event < 256) and **forces CNS off regardless of
`--cns`** — a global CNS setting is therefore safe across mixed ZS/RAW
reprocessing.

### Pulse finding (rewritten 2026-07-24)

Single unified algorithm (replaces the previous legacy-threshold /
derivative-trigger modes; the derivative trigger silently lost 22–46 % of
genuine hits — slow risers never passed its slope threshold, and smoothing lag
placed seeds before the amplitude-threshold crossing, discarding even large
pulses):

1. **Gate** — contiguous runs of samples above `thresholdSigma` (5) × noise RMS.
2. **Bridge** — sub-threshold gaps ≤ `gapMergeSamples` (2) don't end a pulse.
3. **Pile-up split** — recursively split a run at the most prominent interior
   valley (depth ≥ `splitProminenceSigma` (4) × RMS below *both* flanking
   maxima), only if both halves stay ≥ `minSamplesForPeak` wide. A pulse on
   the tail of another makes peak–valley–peak, so pile-up (e.g. nTOF
   gamma-flash regions) is still separated; noise wiggles never reach the
   prominence cut.

Per pulse: local baseline = **median** of up to `baselineLeftWindow` (4)
pre-pulse samples, skipping the `baselineGapSamples` (2) immediately before the
pulse start — those sit on the sub-threshold rise of slow pulses (at a pile-up
split boundary: the valley level, i.e. the previous pulse's tail); amplitude
from a 3-point parabola (or, when ≥
`minSatSamples` (2) consecutive samples sit near ADC max, from least-squares
edge-slope extrapolation); timing at `timingPercentMax` (30 %) of the peak on
the leading edge; threshold crossings interpolated (float); **integral over the
full pulse span** — from where the waveform leaves the baseline to where it
returns, clamped by neighbouring pulses — so the sub-threshold rise and tail
are included (the old between-crossings sum kept only ~half the charge of a
5–10 σ pulse).

### `hits` tree branches

`eventId, trigger_timestamp_ns, channel, amplitude, time, time_of_max, sample,
max_sample, local_baseline, local_max, left_sample, right_sample,
time_over_threshold, integral, saturated, trunc_left, trunc_right,
significance` (+ `feu` after combining). Times are ftst-corrected ns; `left/right_sample` and
`time_over_threshold` are interpolated (not integer samples). `trunc_left` =
pulse already above threshold at the first sample (rise unrecorded — treat
`time` as unreliable); `trunc_right` = still above threshold at the last sample
(tail clipped — `integral`/`amplitude` are floor estimates). The combiner
tolerates older hits files without the trunc branches (defaults false).

### Low-gain / low-threshold mode (nTOF micro-TPC)

For runs at very low gain, `--thr <sigma>` lowers the gate threshold (default
5.0) and `--mf <samples>` enables a **matched-filter gate**: pulse finding
(gate, gap-merge, valley split, crossings/TOT) runs on a boxcar-smoothed copy
of the waveform against the *smoothed*-noise sigma, measured per channel in the
same pedestal pass (`rms_gate` branch of the pedestals tree). Amplitude,
baseline, timing and integral are still measured on the raw waveform.
Rationale: real pulses are many samples wide while the noise tail is dominated
by narrow excursions, so at a fixed fake rate the smoothed gate recovers far
more small pulses than simply lowering the raw threshold — measured on June
det3 data by injecting the measured pulse shape into real pedestal noise:
3-sigma pulses are detected 80 % vs 54 % (raw) at 3 % fakes/waveform, 42 % vs
10 % at 1 % (4-sigma). The noise is time-correlated, so the smoothing gain is
NOT 1/sqrt(k) (boxcar-5 sigma ratio ~0.8); zero-mean/derivative kernels were
tested and lose (the pulse and the correlated noise overlap spectrally).

Recommended operating point: `--thr 3 --mf 5` at 60 ns/sample (`--mf 15` at
20 ns; match the width to the pulse rise duration). Cost: ~2 % fake waveforms
per channel (≈10 isolated fake strips/event on a 512-channel FEU) — downstream
clustering/time-coincidence must absorb that, and on normal-gain data the
added sub-5-sigma hits are noise-dominated (only worth it at low gain). Every
hit carries a `significance` branch (gate peak / gate noise sigma), so a single
low-threshold processing can serve strict analyses too: cutting
`significance >= 5` offline recovers the high-rejection hit set. Default flags
(`--thr 5`, no `--mf`) are byte-identical to the standard analysis.

Investigated and deliberately left unchanged (June det3 data, 2026-07-24, so
this analysis isn't redone): the CNS `nth_element` upper-middle median
convention (offset −0.1 ADC, negligible); MAD/robust-noise thresholds (the
pedestal noise tails are genuine — P(>4·MAD) is ~275× Gaussian — so a
lower robust threshold would roughly double the fake-waveform rate; the
std-based RMS stays); pre-timing smoothing (only ~5 % jitter improvement in a
shape+real-noise MC); the 3-point parabola amplitude (bias ≤2.5 %).

TODO: expose the analyzer configuration (thresholds, split prominence, timing
fraction, …) through the yaml `ConfigManager` in `common/` instead of
compile-time constants in `WaveformAnalyzer.h`.
