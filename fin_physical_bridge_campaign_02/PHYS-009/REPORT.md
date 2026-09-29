# PHYS-009 — VRSPAD-Q12-CALIBRATION-DESK-GATE

Date: 2026-09-29  
Status: **NUMERICAL_EVIDENCE + DESIGN_ONLY**  
Verdict: **PASS_RATE_ENVELOPE_FEASIBILITY__Q12_HARDWARE_VALIDATION_OPEN**

## Frozen target
N=2, label 0 externally anchored, g=3:
`P(j|0)=softmax[(g/2)A7[0,j]]`.

## Public-data desk test
For each public 16-channel VRSPAD rate table, choose 12 distinct physical channels, one control code per channel and one common positive rate scale. Normalize the measured rates and compare with the frozen 12-bin target. No FIN hardware validation record is used.

| operating table | best TV | gate |
|---|---:|---|
| 29V_illum0 | 0.000746463 | PASS |
| 29V_illum2 | 0.000468208 | PASS |
| 29V_illum4 | 0.000739728 | PASS |
| 29V_illum6 | 0.000721472 | PASS |
| 29V_illum8 | 0.000433576 | PASS |
| 29V_illum10 | 0.000326901 | PASS |
| 29V_illum12 | 0.000465636 | PASS |
| 29V_illum14 | 0.000371401 | PASS |
| 29V_illum16 | 0.000634921 | PASS |
| 31V_illum0 | 0.000256719 | PASS |
| 31V_illum1 | 0.000227671 | PASS |
| 31V_illum2 | 0.000238370 | PASS |
| 31V_illum3 | 0.000279710 | PASS |
| 31V_illum4 | 0.000130070 | PASS |
| 31V_illum5 | 0.000255012 | PASS |
| 31V_illum6 | 0.000218990 | PASS |
| 31V_illum7 | 0.000291146 | PASS |
| 31V_illum8 | 0.009152517 | FAIL |

## Main result
- 17/18 operating tables give TV < 0.001.
- One saturated setting, `31V_illum8`, fails at TV≈0.00915 and is retained as a negative operating-point control.
- Selected table: `31V_illum4`, FIN programming TV≈1.3007e-4.
- On that same table all frozen primary countermodel targets are representable with programming TV<=6.12e-4.
- Table-level discretization shifts the FIN-vs-primary-model TV contrasts by at most ≈2.75e-4.

The first FIN bridge therefore only needs one 12-state stochastic unit because the other N=2 label is anchored.

## Remaining hardware gate
This is not a q=12 FIN hardware run. The actual chosen apparatus must still certify q=12 state enforcement, direct 12-bin calibration, drift/cross-talk, invalid states, reset/anchor reproducibility, readout and effective sample size.

Architecture evidence from the public 144-channel code and rate-envelope evidence from the public 16-channel measurements are kept explicitly separate.
