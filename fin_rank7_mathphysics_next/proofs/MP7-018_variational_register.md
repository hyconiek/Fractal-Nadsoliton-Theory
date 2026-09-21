# MP7-018 — variational theorem register after MP7-016/017

Scientific state: **DONE_VARIATIONAL_CLASSIFICATION_FOR_FIRST_TRANSITION**.

The proved variational statements are now separated as follows.

| Statement | Status | Scope |
|---|---|---|
| Exact phase alignment | PROVED_ANALYTIC | Every full-X7 global minimizer has a D12-equivalent aligned nonnegative C4 representative. |
| Boundary support classification | PROVED_ANALYTIC | Exact admissible aligned supports. |
| Nonzero boundary exclusion | PROVED_ANALYTIC | No nonzero aligned boundary stationary root for `0<g<=250/67`. |
| Stationary exhaustion at `g=37/10` | PROVED_INTERVAL_ASSISTED | Exactly uniform + one saddle + one localized minimum. |
| Global minimum at `g=37/10` | PROVED_INTERVAL_ASSISTED | Uniform is unique full-X7 global minimizer. |
| Stationary exhaustion at the equal-energy event | PROVED_INTERVAL_ASSISTED | Only uniform, localized and saddle branch tubes intersect the event parameter box. |
| First global coexistence | PROVED_INTERVAL_ASSISTED | Uniform is unique below `g_eq`; at `g_eq` uniform + 12 localized D12 images minimize globally. |
| Localized branch transversality | PROVED_INTERVAL_ASSISTED | Localized-minus-uniform energy slope is strictly negative at `g_eq`. |
| Orbit selection | NOT PROVED / not supplied | D12 symmetry does not choose one of the 12 localized labels. |

The certified event interval is

`g_eq in [3.7183448971203875, 3.7183448991203876]`.

This closes the bounded global thermodynamic question posed in Package C: the
previous local equal-energy event is now identified as the first global
transition of the supplied finite model.  It does not supply a physical
interpretation of g or select a particular orbit member.
