# 319 — FROZEN-N11-PREDICTOR

Date: 2026-09-27

Status: predictor frozen before any microscopic N=11 construction.

## Goal

Freeze the full observed-process prediction for the next unseen size N=11 before constructing its microscopic state space.

## Preparation of training data

The already-open N=10 ambiguity layer was exactified first. The 82,572 conservatively unlabelled N=10 states reduce to 3,602 D12 orbits. Full deterministic descent of all representatives leaves only 4,452 genuinely unresolved states with stationary mass

`0.0003242666704779515`.

No N=11 information enters this step.

## Frozen N=11 prediction

Clock law: task-318 barrier-aware law

`beta_N = B4 - c/N - d/N^2`

with the task-318 frozen parameters gives

`rho_11(pred) = 0.0011255165589794179`.

The response-shape extrapolation gives

`R_11(pred) = [1.46867369774, 1.73720472044, 0.912327099591, 1, 1.29524474478, 1.27453101879]`.

The second observation is therefore scheduled at

`Delta t = 0.5/rho_11 = 444.2404654209476`.

The preparation map is refit only on already-open N=3..10 data. Any raw prior is projected onto the probability simplex before use. The task-316 rank-1 correction is retained only as an optional diagnostic.

The raw unprojected N=11 preparation forecast again crosses the simplex boundary:

`min p_raw = -0.00920210436`.

Thus simplex sanitation is mandatory.

## Pre-certification against the fixed comparator

Minimum predicted two-time joint-law separation:

`TV(FIN_pred, comparator) = 0.08493141529`.

Frozen task-313 process envelope:

`B_313 = 0.05582267457`.

Therefore the N=11 test is pre-certified before microscopic opening by

`margin = 0.02910874072`.

## Freeze hash

`742d29eb988637df8b170a9cd2a4a2820e3d61a2ad9b1cd7aaeecd2e5f15e98b`

This file and `FROZEN_N11_PREDICTOR_319.json` define the contract used by task 320.
