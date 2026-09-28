# PHYS-013 — TWO-Q12-QUANTIZED-BRIDGE

Date: 2026-09-29

Status: **PUBLIC-CALIBRATION-DRIVEN NUMERICAL DESIGN**

Verdict: **PASS_PREVALIDATION_NUMERICAL_DESIGN__REAL_Q12_DATA_OPEN**

## Construction

Use the public VRSPAD-144 calibration table
`VRSPAD_rates_31.5V_100.0uA.npy` and the public calibration convention
`window_range=0.6 V`, `window_offset=0.07 V`, `center_code=10`.

Select 24 channels using calibration-only criteria:

    |T - 27.25| <= 1.5
    |baseline residual| <= 1%

The public table contains 26 channels satisfying these bounds, so 24 can be
frozen before validation.

Split them into two q=12 nodes.

For a destination channel a with independently calibrated temperature T_a,
program the directed raw weight row

    w_ab = round[T_a E_ab],
    E_ab = -(g/2) A_ab,
    g = 3.

Although the raw integer matrix is slightly asymmetric, the intended
dimensionless interaction `w_ab/T_a` is symmetric. This is temperature
compensation, not a new FIN parameter. The FPGA stores weight rows independently;
the high-level helper's symmetry convention is not a hardware-memory
restriction.

## Numerical result

Using 24 public measured channel calibrations:

- FIN pair-difference TV to the ideal frozen fingerprint:
  **4.13e-4**;
- in 40 random partitions/permutations of the same 24 channels:
  **4.08e-4 to 4.26e-4**;
- raw signed weight codes:
  **-55 to +30**, safely inside int8;
- residual stationary entropy-production / total pair flux:
  about **2.0e-5**, showing that calibration/rounding only weakly breaks
  reversibility in this model.

All four primary PHYS-007 countermodels remain individually representable
inside the 0.003 design envelope under the same temperature-compensation rule.
The closest primary FIN-vs-FLAT contrast changes from 0.011787 ideal to about
0.011612 in the calibrated hardware model.

## Interpretation

This is substantially stronger than a rate-range check. It combines:

- measured per-channel calibration;
- 8-bit weight quantization;
- channel-to-channel temperature variation;
- baseline-rate residuals;
- the full 144-state two-node CTMC;
- the final 12-bin symmetry observable.

It is still **not laboratory validation** because the public data were not
collected from this q=12 / A7 configuration. Cross-talk, simultaneous pulses,
latch invalid states, drift, SPI timing and real closed-loop DAC behavior remain
unmeasured for this configuration.
