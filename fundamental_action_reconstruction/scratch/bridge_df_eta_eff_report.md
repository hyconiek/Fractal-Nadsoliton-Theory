# Scratch bridge audit: D_f-1 vs eta_eff

Status: strong candidate trace; not a legacy→strict bridge theorem.

- `eta_eff` from `bridge_missing_coupling_report.json`: `1.738886298205956`.
- `D_f = 4 ln 2`: `2.772588722239781`; candidate `D_f-1`: `1.772588722239781`.
- Strict `eta`: `1.800000000000000`; midpoint `(eta_eff+eta)/2`: `1.769443149102978`.
- Absolute gaps: `|eta_eff-(D_f-1)|=0.033702424033826`, `|eta-(D_f-1)|=0.027411277760219`, `|midpoint-(D_f-1)|=0.003145573136803`.
- Bracketing: `eta_eff < D_f-1 < eta` is `True`.
- Grep found legacy/path-growth provenance for `D_f=4 ln 2` and `N(d) ~ d^(D_f-1) ≈ d^1.77`; see JSON `grep_audit`.
- No false pass: no kernel identity, no physical-role transfer, no QW-2191 discharge, no ToE closure.
