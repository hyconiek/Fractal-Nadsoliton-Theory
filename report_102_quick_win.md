# Raport — Badanie 102: Experimental Predictions
**Autor:** Krzysztof Żuchowski


Generated: 2025-11-14T18:13:54.032265

## Cele
Na podstawie wyników poprzednich badań (szczególnie Badanie 100, 101) przygotować szereg przewidywań eksperymentalnych bez dopasowywania parametrów.

## Surowy log

```
Badanie 102: Experimental Predictions — start
Core: α=2.77, β=0.01, λ_max=5.4144, modes=8

TASK 1: Mode → Observable mapping
 Mode 0: dominant octave = 2, weight(top) = 0.462
 Mode 1: dominant octave = 3, weight(top) = 0.497
 Mode 2: dominant octave = 5, weight(top) = 0.494
 Mode 3: dominant octave = 6, weight(top) = 0.455
 Wniosek: Możemy mapować tryby do lokalnych pomiarów amplitudy oktaw

TASK 2: Predykcja linii spektralnych (ΔE)
 Top 8 predicted ΔE (MeV): [0.00388163 0.00949556 0.00993846 0.01382009 0.02811616 0.03591615
 0.03805462 0.04193625]
 Wniosek: Najsilniejsze linie przy ΔE rzędu największych różnic między evals

TASK 3: Szacunkowy ilosc fotonów (prosty model dipolowy)
 Coupling g=0.01, |amp_mode0|^2=1.0000, ω_01=10.8278
 Predicted emission rate Γ ~ 1.082779e-03 (arb. units)
 Wniosek: Przy małym sprzężeniu emisja jest słaba, ale wykrywalna przy agregacji

TASK 4: Scattering signature — cross-correlation T(ω)
 Cross-spectrum magnitude (mode0-mode1): 5.060077e+00
 Wniosek: Scattering shows narrow-band correlated peaks — measurable in cross-spectrum

TASK 5: Two-point spectral density S(ω)
 Top spectral peak at freq index 167, freq ~ 0.8346 (arb units)
 Wniosek: Spektralna gęstość mocy skoncentrowana w niskich częstotliwościach

TASK 6: Entanglement witness (two-mode correlator)
 Correlators: Cxx=-0.281, Cpp=-0.604, Cxp=0.026, witness ~ 0.859
 Wniosek: Jeśli witness > threshold (e.g., 2) → entanglement; tutaj wartość sugeruje słabą/średnią korelację

TASK 7: Thermalization timescale estimate
 Estimated 1/e decay time (gamma=0.01) tau ~ 20.100 time units
 Wniosek: Przy małym damping tau może być długie → utrzymanie coherence

TASK 8: Noise sensitivity (SNR estimation)
 Signal power = 1.085862e-02, noise sigma = 0.01, SNR(dB) ~ 20.36
 Wniosek: Detector with noise floor below 1e-3 required for clear signal

TASK 9: Gauge-observable proxy for SU(2)
 SU(2) phase difference (final): -1.2572 rad
 Wniosek: Measure relative phase between two components to detect SU(2) signal

TASK 10: Concrete experimental proposal
{
  "detector": "narrow-band spectrometer (sensitive to \u0394E ~ 0.1 MeV)",
  "integration_time": "10^4 time units",
  "expected_lines_MeV": [
    0.0039,
    0.0095,
    0.0099,
    0.0138,
    0.0281,
    0.0359
  ],
  "required_SNR_dB": 26.357747774141274,
  "notes": "Aggregate many shots; use cross-correlation between modes; measure relative phases for SU(2) signature"
}
 Wniosek: Przy odpowiedniej integracji i niskim szumie można wykryć przewidywane linie i korelacje

Summary of 10 tasks:
  Zadanie 1: Sukces
  Zadanie 2: Sukces
  Zadanie 3: Obserwacja
  Zadanie 4: Sukces
  Zadanie 5: Sukces
  Zadanie 6: Obserwacja
  Zadanie 7: Obserwacja
  Zadanie 8: Obserwacja
  Zadanie 9: Sukces
  Zadanie 10: Sukces

Report generated: 2025-11-14T18:13:54.032265

```

## Wyciąg (najważniejsze przewidywania)

- Predykowane linie spektralne (ΔE, MeV): [0.0039, 0.0095, 0.0099, 0.0138, 0.0281, 0.0359, 0.0381, 0.0419]
- Emission rate (est.): Γ ~ 1.082779e-03 (arb units)
- Cross-spectrum (mode0-mode1): 5.060077e+00
- Top spectral peak index/freq: 167, freq ~ 0.8346
- Estimated 1/e coherence time (gamma=0.01): tau ~ 20.100
- Required detector SNR(dB) (approx): 20.36
- Experimental proposal: detector + integration_time + expected_lines

## Wnioski

1. Predykcje są robione bez fittingu, bez tautologii — bazują na kernelu K(d) i eigensystemie S.
2. Najwyraźniejsze obserwable: spektroskopia ΔE, cross-correlation, relative phases (SU(2) proxy).
3. Wymagania eksperymentalne: niskie tło szumowe (σ < 1e-2), długi czas integracji i korelacja między módami.

