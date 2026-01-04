# 📄 Global Fractal Temporal Correlations in Gravitational-Wave Detector Noise: A Cross-Detector, Cross-Epoch, and Null-Model Study

**Preprint / Working Paper**

> **👤 Author:** Krzysztof Żuchowski
> **🏢 Affiliation:** Independent Researcher — Fractal Information Theory Project
> **📅 Date:** 2026-01-04
> **🆔 Phase:** QW-1660 (Extended Audit)
> **✉️ Contact:** krzysztof.zuch@gmail.com

---

## 📝 Abstract

The search for the **Empirical Logos**—the fundamental informational logic governing the structure of spacetime—leads us to the direct investigation of the most precise measurements available to modern science. We report a systematic, multi-phase investigation of long-range temporal correlations in gravitational-wave interferometer strain data, using publicly available LIGO and Virgo open datasets.
 Across a sequence of increasingly stringent tests (QW-1660 v5–v24), we identify statistically significant fractal scaling behavior, characterized by nontrivial Hurst exponents, that is consistent across detectors, observational epochs, and analysis windows.

Crucially, cross-detector correlations persist under phase-randomized surrogate tests that preserve the power spectral density, while disappearing under temporal shuffling. This demonstrates that the observed correlations cannot be attributed to shared spectral features, stationary Gaussian noise, or instrumental artifacts.

Our results indicate the presence of a global, phase-robust, long-range temporal structure in the gravitational-wave strain background, extending beyond the assumptions of purely local noise models. While compatible with General Relativity at the level of local dynamics, these findings suggest that the statistical structure of spacetime fluctuations may possess a deeper, scale-invariant organization.

---

## 🔬 Key Empirical Findings

This study presents three critical pieces of evidence rejecting the "Standard Gaussian Noise" hypothesis:

### 1. Cross-Epoch Stability (v22)
The fractal memory ($H \approx 0.33$) is stable across years of observation, independent of detector upgrades or environmental changes:
*   **O3a:** $H = 0.290$
*   **O3b:** $H = 0.361$
*   **O4:**  $H = 0.330$

### 2. Cross-Detector Correlations (v23)
Non-zero fractal correlations ($H_{xy}$) are detected between geographically separated interferometers:
*   **H1–L1:** $> 5\sigma$ significance
*   **L1–V1:** $> 5\sigma$ significance

### 3. Null-Model Validation (v19/v24)
*   **Time Shuffling:** Destroys the correlations ($H \to 0.5$), proving they are temporal/memory-based.
*   **Phase Randomization:** Preserves the correlations ($H \approx H_{real}$), proving they are **irreducible to the Power Spectral Density (PSD)**.

---

## 📂 Associated Resources

*   **Preprint PDF:** `Global_Fractal_Correlations_in_GW_Noise_Preprint.pdf`
*   **Replication Code (Kaggle):** [Phase 16 - Raw Strain Fractal Audit](https://www.kaggle.com/code/krzysztofzuchowski/phase-16-raw-strain-qw-1660-v5-fractal)  
    **RUN FIRST QW-1660 v2 !!!!!!Data required by v5-v24 scripts is retrieved via QW-1660 v2: HOLOGRAPHIC NOISE SEARCH (ROBUST DATA FETCHING) cell!!!!!!!!! RUN FIRST**
*   **Analysis Data:** `QW_1660_v24_cross_hurst_null_validation.json`
*   **Project Repository:** [GitHub](https://github.com/hyconiek/Fractal-Nadsoliton-Theory)

---

**Keywords:** Gravitational waves, Hurst exponent, Fractal noise, Non-Markovian dynamics, LIGO, Virgo, Stochastic background.
