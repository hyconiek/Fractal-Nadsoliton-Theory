# 321 — FROZEN-N12-PREDICTOR

Date: 2026-09-27

Status: predictor frozen before any microscopic N=12 construction.

## Clock update without changing architecture

The task-318 architecture is kept fixed:

`beta_N = B4 - c/N - d/N^2`,

with `B4 = 0.66221913712746`.

Only c,d are refit after opening N=11. Using exact rho for N=3..11 gives

`c = 0.7767912392556902`

`d = -2.026973952936632`

and the frozen N=12 forecast

`beta_12(pred) = 0.6115627418626569`

`rho_12(pred) = 0.0005965403181301583`.

## Response shape

Linear response-shape law fitted only to opened N<=11 gives

`R_12(pred) = [1.49297989769, 1.78773629061, 0.907724561791, 1, 1.32281702031, 1.29989534085]`.

The preparation map is deliberately *not* refit using N=11 ambiguous basin labels. Task-319 preparation coefficients trained through exactified N=10 are extrapolated unchanged to N=12, followed by mandatory simplex projection.

Raw minimum prior probability before projection:

`-0.01068332422`.

## Pre-certification

Minimum predicted separation from the fixed comparator:

`0.09039884715`.

Against frozen task-313 envelope

`0.05582267457`

the blind margin is

`0.03457617258`.

Freeze hash:

`248a42aaf5569f5ba50fdeccc380e7b5d13c1220a5c28c04dd4a6c2634eed8c8`.
