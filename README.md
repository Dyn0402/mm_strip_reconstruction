# mm_strip_reconstruction
Reconstruction algorithms build specifically for Micromegas 2D strip detectors built at CEA Saclay.

## Pipeline

`decode` (raw fdf → `nt` sample tree) → `analyze_waveforms` (→ per-FEU `hits`
tree) → `combine_feus_hits` (→ combined hits with `feu` branch) →
`clusterize1d`. `orchestrator/process_run.py` drives the chain for a run.

## Waveform analysis (`waveform_analysis/`)

```
analyze_waveforms <input.root> [output.root] [pedestal.root] [--tps <ns>] [--cns <0|1>]
```

Per event: pedestal subtraction → densify zero-suppressed samples → optional
common-noise subtraction (median per 64-channel block per sample; ON by
default) → pulse finding → one `hits` entry per pulse.

Pedestals: mean = raw baseline per channel; RMS = noise measured **after** CNS
(the data path is CNS-subtracted, so thresholds must be calibrated on post-CNS
noise). Without a pedestal file the data is assumed zero-suppressed
(baseline 256, RMS 1). CNS on zero-suppressed (sparse) data is signal-biased —
the analyzer warns once; prefer `--cns 0` there.

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
pre-pulse samples (at a pile-up split boundary: the valley level, i.e. the
previous pulse's tail); amplitude from a 3-point parabola (or, when ≥
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
time_over_threshold, integral, saturated, trunc_left, trunc_right`
(+ `feu` after combining). Times are ftst-corrected ns; `left/right_sample` and
`time_over_threshold` are interpolated (not integer samples). `trunc_left` =
pulse already above threshold at the first sample (rise unrecorded — treat
`time` as unreliable); `trunc_right` = still above threshold at the last sample
(tail clipped — `integral`/`amplitude` are floor estimates). The combiner
tolerates older hits files without the trunc branches (defaults false).

TODO: expose the analyzer configuration (thresholds, split prominence, timing
fraction, …) through the yaml `ConfigManager` in `common/` instead of
compile-time constants in `WaveformAnalyzer.h`.
