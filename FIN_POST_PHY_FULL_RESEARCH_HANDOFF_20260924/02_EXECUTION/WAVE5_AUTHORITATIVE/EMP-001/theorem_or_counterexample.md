# EMP-001 readiness gate note

The readiness predicate is conjunctive. Let `R_platform` contain the independent calibration, reset, detector-memory, complete-outcome, drift, cost and custody records required by OP-005. EMP-001 may return READY only if every required field has a valid value satisfying the protocol inequalities and the total certified information budget fits the declared platform budget.

Here `R_platform` is absent. Therefore the logical result is `NOT_READY` without assigning synthetic values. Replacing the missing apparatus record by simulation would violate the task's empirical gate.
