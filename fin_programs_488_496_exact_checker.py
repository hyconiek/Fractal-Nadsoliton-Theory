#!/usr/bin/env python3
"""Dependency-minimal exact checks for FIN Programs P488--P496.

The checker uses only the Python standard library.  It replays finite rational
instances of four algebraic identities used by the analytical report:

* the regular weighted-Laplacian Dirichlet identity;
* max-W/min-A Rayleigh complementarity for A=sI-W;
* the block Green/Schur inverse identity;
* the eta=1 damping automaton underlying the Legacy* hyperbolic envelope.

Passing these checks is not a proof-assistant theorem and does not certify the
transcendental frozen strict weights.  It is a dependency-minimal exact replay
of the algebraic core.
"""

from __future__ import annotations

import json
import random
from fractions import Fraction
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "FIN_Programs_488_496_Exact_Check.json"


def zeros(rows: int, cols: int) -> list[list[Fraction]]:
    return [[Fraction(0) for _ in range(cols)] for _ in range(rows)]


def identity(n: int) -> list[list[Fraction]]:
    out = zeros(n, n)
    for i in range(n):
        out[i][i] = Fraction(1)
    return out


def transpose(a: list[list[Fraction]]) -> list[list[Fraction]]:
    return [list(row) for row in zip(*a)]


def matmul(
    a: list[list[Fraction]], b: list[list[Fraction]]
) -> list[list[Fraction]]:
    bt = transpose(b)
    return [
        [sum((x * y for x, y in zip(row, col)), Fraction(0)) for col in bt]
        for row in a
    ]


def matsub(
    a: list[list[Fraction]], b: list[list[Fraction]]
) -> list[list[Fraction]]:
    return [[x - y for x, y in zip(ra, rb)] for ra, rb in zip(a, b)]


def inverse(a: list[list[Fraction]]) -> list[list[Fraction]]:
    n = len(a)
    aug = [row[:] + eye for row, eye in zip(a, identity(n))]
    for col in range(n):
        pivot = next((r for r in range(col, n) if aug[r][col] != 0), None)
        if pivot is None:
            raise ValueError("singular exact matrix")
        aug[col], aug[pivot] = aug[pivot], aug[col]
        scale = aug[col][col]
        aug[col] = [x / scale for x in aug[col]]
        for row in range(n):
            if row == col:
                continue
            factor = aug[row][col]
            if factor:
                aug[row] = [x - factor * y for x, y in zip(aug[row], aug[col])]
    return [row[n:] for row in aug]


def quadratic(a: list[list[Fraction]], x: list[Fraction]) -> Fraction:
    return sum(
        (x[i] * a[i][j] * x[j] for i in range(len(x)) for j in range(len(x))),
        Fraction(0),
    )


def submatrix(
    a: list[list[Fraction]], rows: Iterable[int], cols: Iterable[int]
) -> list[list[Fraction]]:
    rr, cc = list(rows), list(cols)
    return [[a[i][j] for j in cc] for i in rr]


def circulant_regular_weights(
    n: int, radial: list[Fraction]
) -> list[list[Fraction]]:
    w = zeros(n, n)
    for i in range(n):
        for j in range(n):
            if i != j:
                d = min((i - j) % n, (j - i) % n)
                w[i][j] = radial[d - 1]
    return w


def exact_dirichlet_and_resonance_trials(seed: int = 20260809) -> dict:
    rng = random.Random(seed)
    trials = 80
    for _ in range(trials):
        radial = [Fraction(rng.randint(1, 9), rng.randint(2, 13)) for _ in range(3)]
        w = circulant_regular_weights(6, radial)
        s = sum(w[0], Fraction(0))
        a = [[(s if i == j else Fraction(0)) - w[i][j] for j in range(6)] for i in range(6)]
        x = [Fraction(rng.randint(-7, 7), rng.randint(1, 9)) for _ in range(6)]
        lhs = quadratic(a, x)
        rhs = sum(
            (
                w[i][j] * (x[i] - x[j]) ** 2 / 2
                for i in range(6)
                for j in range(6)
            ),
            Fraction(0),
        )
        if lhs != rhs:
            raise AssertionError("Dirichlet identity failed")
        norm = sum((v * v for v in x), Fraction(0))
        if norm:
            rw = quadratic(w, x) / norm
            ra = quadratic(a, x) / norm
            if rw + ra != s:
                raise AssertionError("Rayleigh complementarity failed")
    return {
        "trials": trials,
        "dirichlet_identity": True,
        "rayleigh_complementarity": True,
    }


def exact_schur_trials(seed: int = 20260810) -> dict:
    rng = random.Random(seed)
    trials = 60
    for _ in range(trials):
        # M=R^T R + I is strictly positive and exactly invertible over Q.
        r = [
            [Fraction(rng.randint(-3, 3), rng.randint(1, 5)) for _ in range(4)]
            for _ in range(4)
        ]
        m = matmul(transpose(r), r)
        for i in range(4):
            m[i][i] += Fraction(1)
        e, h = [0, 1], [2, 3]
        mee = submatrix(m, e, e)
        meh = submatrix(m, e, h)
        mhe = submatrix(m, h, e)
        mhh = submatrix(m, h, h)
        schur = matsub(mee, matmul(matmul(meh, inverse(mhh)), mhe))
        retained_inverse = submatrix(inverse(m), e, e)
        if inverse(schur) != retained_inverse:
            raise AssertionError("Schur/Green identity failed")
    return {"trials": trials, "green_schur_identity": True}


def exact_legacy_automaton() -> dict:
    beta = Fraction(1, 100)
    c = Fraction(1)
    states = []
    for d in range(25):
        target = Fraction(1, 1) / (Fraction(1, 1) + beta * d)
        if c != target:
            raise AssertionError("eta=1 damping automaton failed")
        states.append({"d": d, "numerator": c.numerator, "denominator": c.denominator})
        # eta=1 has Delta[(d+1)^eta-d^eta]=1.
        c = c / (Fraction(1) + beta * c)
    return {
        "beta": "1/100",
        "steps": len(states),
        "eta_one_automaton": True,
        "last_state": states[-1],
    }


def run_exact_checks(write: bool = True) -> dict:
    result = {
        "program": "P496",
        "status": "accepted_exact_standard_library_replay",
        "scope": "finite rational algebraic identities; not proof-assistant checked",
        "dirichlet_resonance": exact_dirichlet_and_resonance_trials(),
        "schur": exact_schur_trials(),
        "legacy_automaton": exact_legacy_automaton(),
    }
    if write:
        OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


if __name__ == "__main__":
    print(json.dumps(run_exact_checks(write=True), indent=2, sort_keys=True))
