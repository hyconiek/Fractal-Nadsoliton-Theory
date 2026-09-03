#!/usr/bin/env python3
"""Shared finite OA discrimination functions for ST1272--ST1361."""
from __future__ import annotations

import numpy as np
from scipy.linalg import expm

from fin_total_nadsoliton_common import N, STRICT_A


EIG = np.linalg.eigvalsh(STRICT_A)
DELTA_EIG = EIG[:, None] - EIG[None, :]


def classical_return(t: float) -> float:
    return float(np.mean(np.exp(-EIG * t)))


def quantum_return(t: float) -> float:
    return float(abs(np.mean(np.exp(-1j * EIG * t))) ** 2)


def dephased_quantum_return(t: float, gamma: float) -> float:
    terms = np.exp(-1j * DELTA_EIG * t - 0.5 * gamma * DELTA_EIG**2 * t)
    return float(np.real(np.sum(terms) / N**2))


def classical_distribution(t: float, source: int = 0) -> np.ndarray:
    return np.asarray(expm(-t * STRICT_A)[:, source], dtype=float)


def quantum_distribution(t: float, source: int = 0) -> np.ndarray:
    column = expm(-1j * t * STRICT_A)[:, source]
    return np.asarray(abs(column) ** 2, dtype=float)


def return_derivatives(t: float) -> dict[str, float]:
    z = np.exp(-1j * EIG * t)
    a = np.mean(z)
    ap = np.mean(-1j * EIG * z)
    app = np.mean(-(EIG**2) * z)
    q = float(abs(a) ** 2)
    qp = float(2 * np.real(np.conj(a) * ap))
    qpp = float(2 * (abs(ap) ** 2 + np.real(np.conj(a) * app)))
    c = classical_return(t)
    cp = float(np.mean(-EIG * np.exp(-EIG * t)))
    cpp = float(np.mean(EIG**2 * np.exp(-EIG * t)))
    return {"c": c, "cp": cp, "cpp": cpp, "q": q, "qp": qp, "qpp": qpp}

