# REPLAY STATUS

Fresh handoff replay found one provenance checker with a machine/library-sensitive
over-tight scalar tolerance:

`branch_completion_uniform_76_78_check.py`

Original assertion:
`abs(sfold-1.140334204309495) < 2e-13`

Current SciPy result differs by about:
`4.88498e-13`.

This does not change the reported fold value at the quoted scientific precision.
A non-destructive corrected replay copy using `1e-11` tolerance is included.

Corrected replay return code: 0

Corrected replay output:
```
PASS
pure_k4_fold_t 1.0356584794018928
pure_k4_fold_g 4.993057757778401
pure_k4_fold_s 1.1403342043090066
uniform_thresholds {6: 5.123427551398616, 5: 5.2205547969491946, 4: 5.455614632675739, 3: 6.118057519136048}
k5_D4Phi 0.0550374041045971
k6_D4Phi 0.07619189880372565
```

The original script remains unchanged under `02_REPLAY_CHECKS/`.
