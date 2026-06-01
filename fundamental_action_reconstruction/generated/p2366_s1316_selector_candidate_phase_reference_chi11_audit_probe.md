# P2366 S1316: selector candidate phase-origin chi11 audit

Status: `OPEN_PROGRESS_SELECTOR_CANDIDATE_AUDIT_PHASE_ORIGIN_PREMISE_NO_QW2191_DISCHARGE`

## Result

P2366/S1316 selects the phase-origin selector candidate as the strongest current concrete candidate for `chi11_selector_source`: chiral bispectrum fixes orientation and calibrated coprime Fourier phase recovers source.

## Checks

- Rows checked: `24` / expected `24`.
- Unique orientation/source pairs using mode 1: `24`.
- All coprime modes recover all sources: `True`.
- Translation-invariant magnitudes source-blind: `True`.
- Non-coprime modes alias sources: `True`.

## Hard Limits

- No strict-core phase-origin theorem is claimed.
- No beta_tors -> chi11 theorem or legacy role transfer is claimed.
- QW-2191 remains open; no selector closure or ToE closure is claimed.
