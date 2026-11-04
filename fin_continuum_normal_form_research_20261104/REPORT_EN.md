# Continuum-only quadratic reduction with the correct phase projection

CN implements the next LF task: a source-specific quadratic normal-form
identity for the extended even-continuum response. It uses the whole actual
continuous fields and explicitly retains bound-profile forcing and the
finite-rank terms introduced by continuous projection.

The important bookkeeping point is that physical continuous projection
commutes with the full linear dynamics, but not with the free multipliers
used in half-wave coordinates. The correct projector is J diag(P,P) J^-1,
not componentwise P. With this projector the response decomposition is
exact and its continuous initial phase is zero. This does not alter QP's
separate nonzero whole-state initial chart shift.

The coordinate correction and the newly isolated localized cubic history
both meet the O(epsilon^2) integrated trace budget on critical time. The
extended quadratic source becomes explicit localized quadratic operators
plus retained cubic cascades. Their generic bounds do not yet close the
gate: the localized quadratic history and structured cubic history remain.

Next: establish the stronger mixed local/global estimate for the displayed
localized quadratic operators, then control the explicit cubic history.
The stronger estimate is a target, not a consequence claimed from QP's
existing global localization theorem. No full critical damping law or
physical validation is reported. All new content is in English.
