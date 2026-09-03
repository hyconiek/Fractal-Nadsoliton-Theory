#!/usr/bin/env python3
"""Outward interval helpers for the frozen decimal FIN spectrum.

The certificate model treats the listed decimal spectral values and integer
multiplicities as exact inputs.  This is deliberately narrower than an exact
symbolic certificate for the original kernel parameters.
"""
from __future__ import annotations

import math
import mpmath as mp


DECIMAL_EIGENVALUES = [
    "0",
    "0.75412115420707981",
    "1.5770495144276093",
    "1.9614068619764458",
    "2.1995688493332102",
    "2.2986062720790965",
    "2.3421820411462999",
]
MULTIPLICITIES = [1, 2, 2, 2, 2, 2, 1]


def setup(dps: int = 25):
    iv = mp.iv
    iv.dps = dps
    lambdas = [iv.mpf([s, s]) for s in DECIMAL_EIGENVALUES]
    weights = [iv.mpf([str(m), str(m)]) / 12 for m in MULTIPLICITIES]
    return iv, lambdas, weights


def setup_uncertain(radius: str, dps: int = 25, keep_zero_exact: bool = True):
    """Create symmetric eigenvalue intervals around the frozen decimals."""
    iv = mp.iv
    iv.dps = dps
    rad = iv.mpf([radius, radius])
    lambdas = []
    for index, value in enumerate(DECIMAL_EIGENVALUES):
        center = iv.mpf([value, value])
        if index == 0 and keep_zero_exact:
            lambdas.append(center)
        else:
            lambdas.append(center + iv.mpf([-float(radius), float(radius)]))
    weights = [iv.mpf([str(m), str(m)]) / 12 for m in MULTIPLICITIES]
    return iv, lambdas, weights


def endpoints(x) -> tuple[float, float]:
    return float(x.a), float(x.b)


def abs_upper(x) -> float:
    a, b = endpoints(x)
    return max(abs(a), abs(b))


def scalar_interval(iv, a: str | float, b: str | float | None = None):
    if b is None:
        b = a
    return iv.mpf([str(a), str(b)])


def return_bundle(a: str | float, b: str | float | None = None, dps: int = 25,
                  context=None):
    """Return F=Q-C and its first two derivatives on a time interval."""
    iv, ls, ws = setup(dps) if context is None else context
    t = scalar_interval(iv, a, b)
    x = y = c = x1 = y1 = x2 = y2 = cp = cpp = iv.mpf("0")
    for lam, weight in zip(ls, ws):
        co = iv.cos(lam * t)
        si = iv.sin(lam * t)
        ex = iv.exp(-lam * t)
        x += weight * co
        y += weight * si
        c += weight * ex
        x1 += weight * lam * co
        y1 += weight * lam * si
        x2 += weight * lam * lam * co
        y2 += weight * lam * lam * si
        cp -= weight * lam * ex
        cpp += weight * lam * lam * ex
    f = x * x + y * y - c
    fp = 2 * (-x * y1 + y * x1) - cp
    fpp = 2 * (y1 * y1 - x * x2 + x1 * x1 - y * y2) - cpp
    return f, fp, fpp


def classical_return_interval(t: str | float, dps: int = 25, context=None):
    iv, ls, ws = setup(dps) if context is None else context
    tt = scalar_interval(iv, t)
    value = iv.mpf("0")
    for lam, weight in zip(ls, ws):
        value += weight * iv.exp(-lam * tt)
    return value


def dephased_return_interval(t: str | float, alpha_bounds: tuple[str, str],
                               gamma_bounds: tuple[str, str], dps: int = 20,
                               context=None):
    iv, ls, ws = setup(dps) if context is None else context
    alpha = scalar_interval(iv, *alpha_bounds)
    gamma = scalar_interval(iv, *gamma_bounds)
    tt = alpha * scalar_interval(iv, t)
    total = sum(weight * weight for weight in ws)
    for i in range(len(ls)):
        for j in range(i + 1, len(ls)):
            delta = ls[i] - ls[j]
            total += 2 * ws[i] * ws[j] * iv.exp(
                -iv.mpf("0.5") * gamma * delta * delta * tt
            ) * iv.cos(delta * tt)
    return total


def square_lower(interval_value) -> float:
    a, b = endpoints(interval_value)
    if a <= 0 <= b:
        return 0.0
    return min(a * a, b * b)


def square_upper(interval_value) -> float:
    a, b = endpoints(interval_value)
    return max(a * a, b * b)
