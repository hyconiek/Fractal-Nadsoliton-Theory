# PHYS-010 — STATIC-HARDWARE-VALIDATION-PREREG-FREEZE

Date: 2026-09-29  
Status: **DESIGN_ONLY / NOT EXECUTED**

Freeze only:
- one anchored label i=0;
- one q=12 stochastic node;
- g=3 and g=0 negative control;
- raw 12-bin outcome records;
- independent reset/anchor cycles;
- unchanged PHYS-007 primary countermodels.

Calibration may choose channel permutation, control codes, common rate scale and an operating point using calibration-only records. Validation may change none of these.

Before validation require a direct measured q=12 calibration distribution whose **confidence-certified total TV upper bound is <=0.003**, explicit accounting of invalid one-hot/no-hot/multi-hot events, frozen readout mapping, and a frozen ESS rule if serial correlation exists.

For planning only, a Weissman-type distribution-free 95% bound for 12 bins requires about **628,501 independent calibration samples** for a TV radius of 0.003. A sharper valid preregistered confidence construction may replace this.

No hardware acquisition or lab contact is authorized here. Validation remains unopened.
