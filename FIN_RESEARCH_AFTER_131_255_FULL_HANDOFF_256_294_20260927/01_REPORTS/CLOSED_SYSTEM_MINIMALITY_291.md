# CLOSED-SYSTEM-MINIMALITY-291
## Minimal fresh hidden-content entropy conditionally selects the pure SWAP closure alpha=0

Date: 2026-09-27

Status:
exact finite-horizon information lower bound;
selection of alpha=0 is CONDITIONAL on declaring the visible units to be the complete system and minimizing extra hidden content.

Recall the completion family from report 288.

For n visible units and one-unit rate rho:

- local hidden reset rate per unit:

      alpha rho;

- visible SWAP participation rate per unit:

      (1-alpha) rho.

All alpha can reproduce the same initial one-unit Q3 drift in a uniform product environment.

## 1. Condition on the event schedule

Fix:
- which site/event happens;
- which events are classified as hidden resets versus visible swaps.

Suppose there are

    M_reset

hidden reset events in the observation horizon.

Each exact hidden reset must deliver an independent uniform q-state target symbol.

Conditioned on the schedule, the output target string has entropy

    M_reset log q.

In any globally invertible realization, those symbols must be recoverable from the initial hidden content/controller state.

Therefore:

    boxed:
    H(hidden fresh content | schedule)
      >=
    M_reset log q.

Equivalently, at least

    q^(M_reset)

distinguishable hidden content states are required.

## 2. Poisson mean bound

For q=3, n visible units and local reset rate alpha rho, the expected number of reset events in time T is

    E[M_reset]
      =
    alpha n rho T.

Thus:

    boxed:
    E[H_hidden content | schedule]
      >=
    alpha n rho T log_2(3)

bits.

This lower bound grows linearly in:
- alpha;
- volume n;
- time T.

## 3. Pure SWAP endpoint

For

    alpha=0,

events only permute content already present in the declared visible records.

No fresh target-symbol entropy is required.

The stochastic event SCHEDULE still needs a globally reversible account.

This report does not hide that controller cost.

It compares only the additional fresh CONTENT entropy, conditional on the same declared schedule law.

## 4. Examples at rho=1, n=64, T=100

    alpha=0.01:
      >= 101.44 fresh-content bits

    alpha=0.10:
      >= 1014.38 bits

    alpha=0.50:
      >= 5071.88 bits.

As T->infinity, any fixed alpha>0 needs unbounded fresh-content capacity if exact independent hidden reset behavior is to continue forever.

## 5. Conditional selection theorem

Suppose:

1. the n visible record factors are declared to constitute the complete material/information system under study;
2. the event schedule/controller is held fixed across candidate completions;
3. additional hidden fresh-content entropy is minimized;
4. visible content transport by SWAP is allowed.

Then:

    boxed:
    alpha=0

is the unique minimizer of fresh hidden-content demand.

This gives a concrete source principle for the conservative closure.

## 6. Boundary

This is NOT yet an unconditional FIN theorem.

The decisive premise is:

    "the declared visible record factors are the complete closed system."

If nature contains an additional reservoir, alpha>0 is allowed.

The exact natural extension itself can encode an unbounded history reservoir.

Therefore report 291 does not prove the universe must choose pure SWAP.

It proves:

    boxed:
    once complete closure + minimal hidden content are adopted,
    the reset/SWAP allocation ambiguity is removed.

## Importance

The missing composition law has been reduced to a clear ontological/operational declaration:

    what degrees of freedom count as part of the complete closed system?
