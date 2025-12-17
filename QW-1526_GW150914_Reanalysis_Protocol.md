# QW-1526: GW150914 Reanalysis Protocol
**Objective:** Test for non-standard amplitude scaling in Gravitational Wave data.
**Method:** Bayesian Inference with free amplitude exponent parameter ($n$).

---

## 1. Executive Summary
This study re-analyzes the GW150914 event data to determine if the amplitude scaling follows the General Relativity prediction ($A \propto 1/D_L^{1.0}$) or if the data allows/prefers a non-standard scaling ($A \propto 1/D_L^n$) predicted by Active Vacuum theories like FIN ($n \approx 0.66$).

**Key Constraint:** The phase evolution, polarization, and spin dynamics are kept strictly standard (GR). Only the luminosity distance scaling in the amplitude is modified.

---

## 2. Analysis Configuration

### Data Config
* **Event:** GW150914
* **Detectors:** H1, L1
* **Duration:** 4s centered on coalescence
* **Sampling:** 4096 Hz
* **Bandpass:** 20 Hz - 512 Hz
* **PSD:** Estimated locally (Welch method, 4s segments)

### Waveform Model
* **Baseline:** `IMRPhenomD` (Frequency Domain)
* **Modification:**
  $$ \tilde{h}_{FIN}(f) = \frac{1}{D_L^n} \times \mathcal{M}_{IMR}(f, \theta_{GR}) \times e^{i \Psi_{IMR}(f, \theta_{GR})} $$

---

## 3. Priors

| Parameter | Prior Distribution | Bounds |
| :--- | :--- | :--- |
| Chirp Mass $\mathcal{M}_c$ | Uniform (Detector Frame) | [20, 50] $M_\odot$ |
| Mass Ratio $q$ | Uniform | [0.125, 1.0] |
| Luminosity Dist $D_L$ | Power Law ($D_L^2$) | [100, 2000] Mpc |
| **Scaling Exp $n$** | **Uniform** | **[0.5, 1.2]** |
| Spin Magnitudes | Uniform | [0, 0.99] |
| Inclination $\iota$ | Sinusoidal | $[0, \pi]$ |

---

## 4. Run Plan & Criteria

### Run A: Reference (GR)
* Fix $n=1$.
* Calculate $\ln \mathcal{Z}_{GR}$.

### Run B: Test (FIN)
* Free $n \in [0.5, 1.2]$.
* Calculate $\ln \mathcal{Z}_{FIN}$ and $P(n|d)$.

### Interpretation Criteria
* **Null Result:** $n=1$ falls within the 90% High Posterior Density (HPD) interval.
* **Tension:** $n=1$ falls outside 90% HPD.
* **Discovery:** $Bayes Factor (\mathcal{Z}_{FIN} / \mathcal{Z}_{GR}) > 10$.

---

## 5. Software Requirements
* `bilby` >= 1.1.0
* `lalsuite`
* Python 3.8+
* Custom waveform wrapper (provided in `fin_waveform.py`)
