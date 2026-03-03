# RAPORT QW-1724: GW METHOD AUDIT

- Data UTC: 2026-03-02T16:14:17.708016+00:00
- Risk points: 17
- Werdykt: **GW_PIPELINE_METHOD_HIGH_RISK_INCONCLUSIVE**

## Zidentyfikowane problemy
- [1] (critical) Sprzecznosc modelu projekcji: v50 daje Cross_H=0.509 (daleko od 0.31), a v55 twierdzi projekcje do 0.301~0.31.
- [2] (critical) Brak inwariancji skali: Rozrzut H_cross po skalach = 1.023 (od -0.131 do 0.892), mimo claimu o inwariancji.
- [3] (major) Niespojnosc fazowa: Po whitening H_cross=0.686 sugeruje silna strukture, ale test opoznienia ma max korelacje -0.0057 przy 73.24 ms (bardzo slabo).
- [4] (major) Brak reprodukcji miedzy parami detektorow: H1-L1=0.311 vs H1-V1=0.715; roznica 0.404.
- [5] (major) Mozliwy drift protokolu null-testu: W v61 obserwacja jest hard-coded jako 0.311 mimo innej dlugosci sygnalu (N=131072).
- [6] (minor) Silna zaleznosc pasmowa: H_cross po pasmach: min=0.413, max=0.731.
- [7] (major) Niska identyfikowalnosc mechanizmu: Zmiana envelope/fazy przesuwa H_cross o 0.419 (Test1=0.287, Test2=0.705).

## Artefakty
- JSON szczegolowy: `report_qw1724_gw_method_audit.json`
