# FIN RESEARCH PROPOSAL: Pulsar Timing Arrays (PTA) Search

**Date:** 2026-02-27
**Target:** NANOGrav 15-yr Data Release / EPTA / PPTA
**Objective:** Search for the geometric signature of the FIN background ($H_{true} \approx 0.23$) in the long-term phase memory of spacetime, bypassing terrestrial hardware PSD artifacts.

## Rationale
The investigation into raw LIGO data (Phases 47-66) proved that terrestrial interferometers share an overwhelming physical design spectral density envelope (PSD). This PSD "tube" corrupts Cross-DFA estimators, generating false correlated fractalities tied to suspension resonances rather than gravitational fields.

To test FIN properly, we desperately need:
1. **Independent hardware topologies** (pulsars are distinct, natural clocks).
2. **No shared electronic/mechanical control loops**.
3. **Macro-scale astrophysical baselines** tracing deep spacetime.

**Pulsar Timing Arrays (PTAs)** fundamentally fulfill all these requirements. 

## Methodology Overview
1. **Data Acquisition:** Download the public NANOGrav 15-yr data release, specifically the decoupled "timing residuals" across ~60-70 millisecond pulsars.
2. **Cross-MF-DFA Implementation:** 
   - Instead of correlating H1/L1, we correlate `Residuals(Pulsar_A)` vs `Residuals(Pulsar_B)`.
   - The analysis scales will be in the regime of months to decades (nHz), perfectly tracking the slow fractional drift expected in a spatial memory lattice.
3. **Robust Null Testing:**
   - **Phase-randomization:** Does $H \approx 0.23$ survive if we destroy the phases of the residuals?
   - **Hellings-Downs Isolation:** We must distinguish standard GW stochastic background (quadripolar correlation) from the FIN field structure.
   - **Envelope Swap (Red Noise Null):** Pulsars have intrinsic spin-down red noise. We will swap true phases into generic red-noise power laws to ensure we aren't just measuring independent pulsar aging.

## Deliverables
- Python/Jupyter pipeline optimized for non-uniform time series (PTA observations aren't perfectly uniform like LIGO's 4096 Hz ticks, requiring Lomb-Scargle or interpolated DFA adaptations).
- Cross-correlation fractional exponent matrix for all pulsar pairs.

The hunt for $H=0.23$ moves to the stars.
